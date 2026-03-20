# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# function to get the file address using location of files and file name

# get file path; file_name if your file name;
# memade set to TRUE if in a folder you designated as made by you (labeled output_data);
# change data_folder and output_data in the constats file
#'
#' @description
#' Get file path; file_name if your file name.
#'
#' @param file_name name of file to access
#' @param memade set to FALSE. Set to TRUE to if you have a separate folder with datasets you created in the data folder called "memade".
#' @returns Full file path with name is returned.
#' @export
getfp = function(file_name, memade = FALSE) {
  return(
    if(memade == FALSE) {
      paste0(data_folder, "/", file_name)
    } else {paste0(output_data, "/", file_name) }
  )
}

#' @description
#' Function to summarise variables quickly
#' @param dataframe dataframe
#' @param to_summarise variables to summarise
#' @param ... arguments to pass on to group_by.
#' @returns A summary tables with mean, standard deviation, number of observations, and number of missing observations.
#' @export
summarise2 = function(dataframe, to_summarise, ...) {

    df = dataframe %>%
      group_by(!!! ensyms(...)) %>%
      summarise(
        mean = mean({{to_summarise}}, na.rm=TRUE),
        sd = sd({{to_summarise}}, na.rm=TRUE),
        n = sum(!is.na({{to_summarise}})),
        n_na = sum(is.na({{to_summarise}}))
      )

  return(df)
}

#' @description
#' Function to efficiently import datasets
#' @param file_path file_path
#' @param ... columns to importa, except for ID and date columns
#' @param day name of column with number denoting day.
#' @param month name of column with number denoting month.
#' @param year name of column with number denoting year.
#' @param name_date name to give the new full date columnt. Default is set to "date".
#' @returns Imported dataset with specified columns and a column of dates.
#' @export
import.data <- function(file_path,
                        ..., # which columns, besides ID (uncchdid) and date, do you want returned?
                        day = "none", month = "none", year = "none", # date
                        name_date = "date" # need to give date columns names dep on the file name
                        ) {
  # reading in file
  if (grepl("dta$", file_path, ignore.case = TRUE)) {
    df = haven::read_dta(file_path)
  } else {
    df = get(load(file_path))
  }

  # checking if file has uncchdid
  if (!("uncchdid" %in% colnames(df)) & !("uncid" %in% colnames(df))) {
    df = left_join(df, dplyr::select(haven::read_dta(getfp("idbaselinkid.dta")), -outcome), by = c("basebrgy", "basewman")) %>%
      dplyr::select(-basebrgy, -basewman)
  } else if (!("uncchdid" %in% colnames(df))) {
    df = left_join(dplyr::filter(df, uncid >= 20000), dplyr::select(haven::read_dta(getfp("idbaselinkid.dta")), -outcome), by = c("basebrgy", "basewman")) %>%
      dplyr::select(-basebrgy, -basewman, -uncid)
  } else if ("basebrgy" %in% colnames(df) & "basewman" %in% colnames(df)) {
    df = dplyr::select(df, -basebrgy, -basewman)
  }
  # getting date from file
  if (year != "none") {
    if (median(df[[year]], na.rm = TRUE) >= 1900) {
      df = df %>% mutate(
        date = as.Date(paste0(eval(parse(text = day)),"/", eval(parse(text = month)), "/", eval(parse(text = year))), format = '%d/%m/%Y')
      ) %>% rename(!!name_date := date) %>% dplyr::select(uncchdid, !!name_date, ...)
    } else {
      df = df %>% mutate(
        date = as.Date(paste0(eval(parse(text = day)),"/", eval(parse(text = month)), "/", eval(parse(text = year))), format = '%d/%m/%y')
      ) %>% rename(!!name_date := date) %>% dplyr::select(uncchdid, !!name_date, ...)
    }
  } else {
    df = dplyr::select(df, uncchdid, ...)
  }
  return(df)
}

#' @description
#' Residualises epigenetic clocks estimates on chronological age.
#' @param dataframe dataframe
#' @param clock epigenetic clock column name
#' @param age chronological age column name
#' @param outptut optional parameter. Returns linear regression residuals of epigenetic age estimates on age if left empty.
#' @returns If output == plot, returns plot of the linear regression. Else, outputs a dataframe with the newly-made residuals.
#' @export
make.residuals0 = function(dataframe, clock, age, output) {
  if (is.character(clock)) {
    column_name <- clock
    # Check if the column name exists in the dataframe
    if (!(column_name %in% colnames(dataframe))) {
      stop(paste("Column", column_name, "does not exist in the dataframe."))
    }
    # Convert column name to symbol
    clock <- sym(column_name)

  } else {
    # If the input is not a string, assume it's a column object
    column_name <- deparse(substitute(column_input))
    # Check if the column object exists in the dataframe
    if (!(column_name %in% colnames(data))) {
      stop(paste("Column", column_name, "does not exist in the dataframe."))
    }
  }
  dataframe = dataframe %>% filter(!is.na({{clock}})) %>% filter(!is.na({{age}}))
  mod = lm(reformulate(deparse(substitute(age)), deparse(substitute(clock))), data = dataframe)
  dataframe[[paste0("rs.",column_name)]] = rstudent(mod)
  dataframe[[paste0("rst.",column_name)]] = rstandard(mod)
  dataframe[[paste0("r.",column_name)]] = residuals(mod)

  if (output == "plot") {
    return(plot(mod))
  } else if (output == "student") {
    return(dplyr::select(dataframe, uncchdid, {{clock}}, paste0("rs.",column_name), {{age}}))
  } else if (output == "all") {
    return(dplyr::select(dataframe, uncchdid, {{clock}}, paste0("r.",column_name), paste0("rs.",column_name), paste0("rst.",column_name), {{age}}))
  } else if (output == "standard") {
    return(dplyr::select(dataframe, uncchdid, {{clock}}, paste0("rst.",column_name), {{age}}))
  } else {
    return(dplyr::select(dataframe, uncchdid, {{clock}}, paste0("r.",column_name), {{age}}))
  }
}

cor.res = function(dataframe, clock, which = "standard") {
  zscore = paste0("z.r.", clock)
  if (which == "standard") {
    rs = paste0("rs.", clock)
  } else if (which == "student") {
    rs = paste0("rst.", clock)
  } else {
    stop()
  }
  return(cor.test(dataframe[[zscore]], dataframe[[rs]], na.rm = TRUE))
}

#' @description
#' Creates linear models and, optionally, presents them as a table.
#' @param output outcome
#' @param exposure exposure
#' @param covars covariates
#' @param dataframe dataframe
#' @param to_return what to return. Default is set to "table"
#' @param my_type what type of table to return; goes into stargazer.
#' @param ... arguments to be passed on to stargezer
#' @returns If output == "table", returns a regression coefficient table. It output is anything else, returns a list of linear regressions.
#' @export
table.mvreg = function(output, exposure, covars, dataframe, to_return = "table", my_type = "html", ...) {
  mvreg.list = list()
  for (i in 1:length(output)) {
    mvreg.list[[i]] <- lm(reformulate(
      response = output[i],
      termlabels = c(exposure, covars)), data = dataframe
    )
  }

  exposure_name = ifelse(exposure == "z.bmi", "BMI (kg/m2)",
                         ifelse(exposure == "z.waist", "Waist (cm)",
                                ifelse(exposure == "z.whtr", "Waist-to-height ratio", exposure)))
  if ((TRUE %in% grepl(".*\\*sex", covars)) == TRUE) {
    to.keep = c(exposure, "sex", covars[which(grepl(".*\\*sex", covars))])
    covar.labels = c(exposure_name, "Male sex", paste0(sub(" \\(.*", "", exposure_name), "*Male sex"))
  } else if (("sex" %in% covars) == TRUE) {
    to.keep = c(exposure, "sex")
    covar.labels = c(exposure_name, "Male sex")
  } else {
    to.keep = c(exposure)
    covar.labels = c(exposure_name)
  }

  if (to_return == "table") {
    return(
      stargazer(mvreg.list,
                type=my_type,
                header = FALSE,
                digits = 3,
                keep = to.keep,
                style = 'default',
                dep.var.labels = '',
                report="vc*p", star.cutoffs=c(0.05, 0.01, 0.001),
                # ci=T,
                covariate.labels = covar.labels,
                notes = 'signif. levels: * p < 0.05; ** p < 0.01; *** p < 0.001', notes.append = FALSE, model.numbers = FALSE,
                ...,
                digits.extra = 0,
                font.size = "footnotesize", column.sep.width = "3pt",
                omit.stat = "all")
    )
  } else {
    return(mvreg.list)
  }
}

#' @description
#' Creates coefficient plots.
#' @param models list of models whose coefficients to plot
#' @param model_names model names
#' @param to_omit regression coefficients to omit
#' @param legend_position legend position in ggplot2::theme
#' @param xlim if you want to choose x-axis limits, put here a vector (xmin, xmax). Otherwise, chosen automatically
#' @param to.return set to "plots". If set to "x", it will return (xmin, xmax)
#' @param my.title plot title
#' @param x.lab x-axis label
#' @param ... arguments to be passed on to jtools::plot_summs
#' @returns If to.return == "plots, returns a regression coefficient plot. It to.return is "x", returns c(xmin, xmax).
#' @export
plot.coeffs = function(models, model_names, to_omit = c(), legend_position = "none", xlim = "choose", to.return = "plots", my.title = "choose", x.lab = NULL, ...) {

  # models = all.unadj
  # model_names = mvreg_response_labels
  # to_omit = c("(Intercept)", "sex1")

  if (is.numeric(xlim) == FALSE | length(xlim) < 2) {
    conf.int.low = c()
    conf.int.high = c()

    for (i in 1:length(models)) {
      for (j in 1:length(models[[1]])) {
        conf.int.low <- c(conf.int.low, min(confint(models[[i]][[j]])[2,1]))
        conf.int.high <- c(conf.int.high, max(confint(models[[i]][[j]])[2,2]))
      }
    }

    xmin = min(conf.int.low)
    xmax = max(conf.int.high)
  } else {
    xmin = xlim[1]
    xmax = xlim[2]
  }

  if (to.return == "x") {
    return(c(xmin, xmax))
  } else {

    exposure = c()
    for (i in 1:length(models)) {
      exposure = c(exposure, rownames(as.data.frame(models[[i]][[1]]$coefficients))[2])
    }

    plots = list()
    for (i in 1:length(models)) {
      plots[[i]] = jtools::plot_summs(models[[i]],
                                      robust = TRUE,
                                      point.size = 6, model.names = model_names,
                                      colors = c("#004949","#18BC9C","#FF6DB6","#006DDB","#920000","#FF6E00","#B66DFF","#490092"),
                                      !!! ensyms(...),
                                      omit.coefs = to_omit,
                                      coefs = setNames(exposure[i], "")) +
        ggplot2::labs(
          title = ifelse(my.title == "choose", exposure_labels[i], ""),
          x = ifelse(is.null(x.lab), "Estimate", ifelse(x.lab == "", "", x.lab)),
          y = exposure_labels[i]
        ) +
        ggplot2::theme(
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.position = legend_position
        ) +
        geom_vline(xintercept = 0, color = "red3", linetype = "dashed", size = 0.6) +
        ggplot2::xlim(xmin, xmax) +
        ggpubr::font("title", size = 20, face = "bold") +
        ggpubr::font("legend.title", size = 20) +
        ggpubr::font("legend.text", size = 18) +
        ggpubr::font("xlab", size = 15) +
        ggpubr::font("ylab", size = 20, angle = 90,

                     margin = margin(t = 0, r = -3, b = 0, l = 0)) +
        ggpubr::font("xy.text", size = 20, color = "gray30")
    }

    return(plots)
  }

}

#' @description
#' Function to extract coefficients, p-values, and confidence intervals
#' @param model linear regression model
#' @param model_name optional model name
#' @param to_return default set to "everything"; otherwise returns specified regression coefficients
#' @returns A tidy table with prespecified regression coefficients and their statistics
#' @export
extract_coef_pval <- function(model, model_name = "", to_return = "everything") {

  tidy_model = broom::tidy(model) %>%
    dplyr::select(term, estimate, p.value) %>%
    mutate(term = ifelse(grepl("exposure", term), "Exposure", term)) %>%
    # left_join(
    #   confint(model) %>%
    #     as.data.frame() %>%
    #     mutate(term = rownames(.)), by = "term"
    # ) %>%
    mutate(
      term = case_when(
        grepl("uw1", term) ~ gsub("uw1", "uw", term),
        grepl("sex1", term) ~ gsub("sex1", "Sex", term),
        grepl("sex", term) ~ gsub("sex", "Sex", term),
        TRUE ~ term
      )
    ) %>%
    mutate(
      term = case_when(
        grepl("bmi$", term) ~ "BMI",
        grepl("bmi:", term) ~ gsub(".*bmi:", "BMI*", term),

        grepl("waist$", term) ~ "WC",
        grepl("waist:", term) ~ gsub(".*waist:", "WC*", term),

        grepl("whtr$", term) ~ "WHtR",
        grepl("whtr:", term) ~ gsub(".*whtr:", "WHtR*", term),
        TRUE ~ term
      ),
      model = model_name
    ) %>%
    dplyr::select("Clock" = model, term, estimate, p.value)

  if(paste0(to_return, collapse = ".") == "everything") {
    tidy_model2 = tidy_model
  } else {
    tidy_model2 = tidy_model %>% filter(
      term %in% to_return
    )
  }

  return(as.data.frame(tidy_model2))
}

#' @description
#' Adds Benjamini-Hochberg adjustment
#' @param dataframe dataframe
#' @param estimate estimate
#' @param group which term to adjust
#' @param group2 group for which to adjust, such as set of epigenetic clocks used; default set to "all".
#' @param p.val p-value column name; default set to p.value
#' @param new_name optional what to name "Estimate" column
#' @param keep_all set to FALSE and removes original estimate and p-value columns
#' @returns A tidy table prepared to html output with Benjamini-Hochberg statistical significance adjustment
#' @export
adj_benj_hoch = function(dataframe, estimate = estimate, group = term, group2 = Clock, which.group2 = "all", p.val = p.value, new_name = "Est.", keep_all = FALSE) {

  if (paste(which.group2, collapse = ".") == "all") {
    n = nrow(dataframe)/length(unique(dataframe[[deparse(substitute(!!group))]]))
    dataframe2 = dataframe %>% group_by({{group}}) %>%
      mutate(
        !!new_name := case_when(
          {{p.val}} < 0.001*rank/n ~ paste0(sprintf("%.3f", {{estimate}}), "***<br><i>(p-val &#60; 0.001)</i>"),
          {{p.val}} < 0.01*rank/n ~ paste0(sprintf("%.3f", {{estimate}}), "**<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          {{p.val}} < 0.05*rank/n ~ paste0(sprintf("%.3f", {{estimate}}), "*<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          TRUE ~ paste0(sprintf("%.3f", {{estimate}}), "<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>")
        )
      ) %>% dplyr::select(-rank)
  } else if (paste(which.group2, collapse = ".") == "none") { # if you want just the formatting, and no adjustment
    dataframe2 = dataframe %>%
      mutate(
        !!new_name := case_when(
          {{p.val}} < 0.001 ~ paste0(sprintf("%.3f", {{estimate}}), "***<br><i>(p-val &#60; 0.001)</i>"),
          {{p.val}} < 0.01 ~ paste0(sprintf("%.3f", {{estimate}}), "**<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          {{p.val}} < 0.05 ~ paste0(sprintf("%.3f", {{estimate}}), "*<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          TRUE ~ paste0(sprintf("%.3f", {{estimate}}), "<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>")
        )
      )
  } else {
    n = length(which.group2)
    # subsetting
    dataframe2 = rbind(
      dataframe %>% filter({{group2}} %in% which.group2) %>% group_by({{group}}) %>%
      mutate(
        rank = rank({{p.val}})
      ) %>% ungroup() %>%
      mutate(
        !!new_name := case_when(
          {{p.val}} < 0.001*rank/n ~ paste0(sprintf("%.3f", {{estimate}}), "***<br><i>(p-val &#60; 0.001)</i>"),
          {{p.val}} < 0.01*rank/n ~ paste0(sprintf("%.3f", {{estimate}}), "**<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          {{p.val}} < 0.05*rank/n ~ paste0(sprintf("%.3f", {{estimate}}), "*<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          TRUE ~ paste0(sprintf("%.3f", {{estimate}}), "<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>")
        )
      ) %>% dplyr::select(-rank),
    # the rest, not adjusted
    dataframe %>% filter(!({{group2}} %in% which.group2)) %>%
      mutate(
        !!new_name := case_when(
          {{p.val}} < 0.001 ~ paste0(sprintf("%.3f", {{estimate}}), "***<br><i>(p-val &#60; 0.001)</i>"),
          {{p.val}} < 0.01 ~ paste0(sprintf("%.3f", {{estimate}}), "**<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          {{p.val}} < 0.05 ~ paste0(sprintf("%.3f", {{estimate}}), "*<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>"),
          TRUE ~ paste0(sprintf("%.3f", {{estimate}}), "<br><i>(p-val = ", sprintf("%.3f", {{p.val}}), ")</i>")
        )
      )
    )
  }

  if (keep_all == FALSE) {
    dataframe2 = dataframe2 %>% dplyr::select(-{{estimate}}, -{{p.val}})
  }
  return(dataframe2)
}

#' @description
#' Makes html table
#' @param model_list list of models whose regression coefficients to report
#' @param model_names model names
#' @param keep_betas which coefficients to keep
#' @param my_header optional header name
#' @param ... to pass on to kable.
#' @param by_sex are models stratified by sex and to be reported as such?
#' @param kable default TRUE to return kable output html table
#' @param bha which groups to apply Benjamini-Hochberg to. Default is mvreg_response_labels[2:7].
#' @returns A tidy table prepared to html output with Benjamini-Hochberg statistical significance adjustment
#' @export
make_table2 = function(model_list, model_names, keep_betas, my_header = "Betas", ..., by_sex = FALSE, kable = TRUE, bha = mvreg_response_labels[2:7]) {

  col_length = if (by_sex == FALSE) {
    length(model_list)
  } else {
    length(model_list)/2 # make sure to first put all female models and then all male models
  }

  if (by_sex == FALSE) {
    new_col_name = c("Model")
    new_num = -1
  } else {
    new_col_name = c("F Model", "M Model")
    new_num = 0
  }

  my_list = list()

  for (i in 1:col_length) {
    assign(paste0("mod", i), extract_coef_pval(model_list[[i]][[1]], model_names[1], keep_betas))

    for (j in 2:length(model_names)) {
      assign(paste0("mod", i), rbind(
        get(paste0("mod", i)),
        extract_coef_pval(model_list[[i]][[j]], model_names[j], keep_betas)
      ))
    }

    my_list[[i]] = adj_benj_hoch(get(paste0("mod", i)), new_name = paste(new_col_name[1], new_num + i), which.group2 = bha)

  }

  if (by_sex == TRUE) {
    for (i in 1:col_length) {
      assign(paste0("mod", i + col_length), extract_coef_pval(model_list[[i + col_length]][[1]], model_names[1], keep_betas))

      for (j in 2:length(model_names)) {
        assign(paste0("mod", i + col_length), rbind(
          get(paste0("mod", i + col_length)),
          extract_coef_pval(model_list[[i + col_length]][[j]], model_names[j], keep_betas)
        ))
      }

      my_list[[i]] = full_join(
        adj_benj_hoch(get(paste0("mod", i)), new_name = paste(new_col_name[1], new_num + i), which.group2 = bha),
        adj_benj_hoch(get(paste0("mod", i + col_length)), new_name = paste(new_col_name[2], new_num + i), which.group2 = bha),
        by = c("Clock", "term")
      )

    }
  }

  # combining tables in the list above
  all.table = Reduce(
    function(x, y) full_join(x, y, by = c("Clock","term"), relationship = "many-to-many"),
    my_list[2:col_length],
    init = my_list[[1]]
  ) %>% arrange(Clock) %>%
    mutate(across(everything(), ~ifelse(is.na(.), "", .)))

  kable1 = kable(all.table, format='html', ..., escape=FALSE) %>%
    kable_styling(bootstrap_options = "striped", full_width = F) %>%
    column_spec(column = 1, bold = TRUE) %>%
    column_spec(column = 3:ncol(all.table), width = "16em") %>%
    collapse_rows(columns = 1, valign = "top") %>% add_header_above(header = c(setNames(object = ncol(all.table), nm = my_header)))

  if (kable == TRUE) {
    return(kable1)
  } else {
    return(all.table)
  }

}

#' @description
#' Makes publication-ready html table with BH adjustment.
#' @param model_list list of models whose regression coefficients to report
#' @param model_names model names
#' @param keep_betas which coefficients to keep
#' @param start start of which to BH adjust
#' @param stop stop of which to BH adjust
#' @param which.group2 defaults to list-of-model_names[start:stop]
#' @returns A tidy table prepared to html output with Benjamini-Hochberg statistical significance adjustment
#' @export
# produce tables
make_table2_pub = function(model_list, model_names, keep_betas, start = 1, stop = 7, which.group2 = model_names[start:stop]) {

  # model_list = all.pc.final,
  # model_names = mvreg_response_pc,
  # keep_betas = c("BMI", "WC", "WHR", "WHtR", "BMI*Sex", "WC*Sex", "WHR*Sex", "WHtR*Sex")
  # col.names = c(
  #   "Measure", "Coeff", mvreg_response_labels_pub[5:10]
  # )

  my_table = do.call(
    rbind,
    lapply(start:stop, function(y)
      do.call(
        rbind,
        lapply(
          1:length(exposures),
          function(x) extract_coef_pval(
            model_list[[x]][[y]],
            model_name = exposures_pub[x],
            to_return = keep_betas
          )
        )
      ) %>%
        mutate(across(c(estimate, p.value), as.numeric)) %>%
        arrange(Clock) %>% mutate(EC = model_names[y])
    )
  ) %>% adj_benj_hoch(group2 = EC, which.group2 = which.group2)

  all.table = Reduce(
    function(x, y) full_join(x, y, by = c("Clock","term")),
    do.call(list, lapply((start+1):stop, function(x) my_table[my_table$EC %in% model_names[x],] %>% rename(!!model_names[x] := Est.) %>% dplyr::select(-EC))),
    init = my_table[my_table$EC == model_names[start],] %>% rename(!!model_names[start] := Est.) %>% dplyr::select(-EC)
  ) %>% arrange(Clock)

  return(all.table)

}

my_kable = function(table, ..., my_header = "Betas") {
  kable(table %>% mutate(across(everything(), ~ifelse(is.na(.), "", .))), format='html', escape=FALSE, ...) %>%
    kable_styling(bootstrap_options = "striped", full_width = F) %>%
    column_spec(column = 1, bold = TRUE) %>%
    column_spec(column = 3:ncol(table), width = "16em") %>%
    collapse_rows(columns = 1, valign = "top") %>% add_header_above(header = c(setNames(object = ncol(table), nm = my_header)))
}


#' Calculates summary stats for Stata-like margins plot
#'
#' @description
#' Calculates mean and 95\%CI using t.test function around a set of x-values.
#' Note: Mistral AI was used to help create this function.
#'
#' @param data Dataset
#' @param x_col independent variable (x-axis)
#' @param y_col dependent variable (y-axis)
#' @param x_value x_col value for which to calculate mean and confidence interval
#' @param bandwidth interval around value for which to calculate summary. For example, if bandwidth = 1, summary is calculated for (x_col-1, x_col+1) interval.
#' @returns a dataframe containing x_value, mean of y_col for chosen interval, and lower and upper 95\%CI bounds.
calculate_stats = function(data, x_col, y_col, x_value, bandwidth) {
  if (is.character(x_col)) {
    x_col = sym(x_col)
  }
  data = data %>% filter(!is.na({{x_col}})) %>% filter(!is.na({{y_col}}))
  if (nrow(data) == 0) {
    stop("No non-NA x or y values. Verify that x and y columns are non-empty and ")
  }
  subset = data %>% filter(abs({{x_col}} - x_value) < bandwidth)
  mean_y = subset %>% dplyr::select({{y_col}}) %>% colMeans(na.rm = TRUE)
  ci <- t.test(subset %>% dplyr::select({{y_col}}))$conf.int
  return(data.frame(my_x = x_value, mean = mean_y, lower = ci[1], upper = ci[2]))
}

#' Produces and plots predictive margins
#'
#' @description
#'
#' Produce predictive margins given a model object, new data, and, optionally, plot them.
#' For plotting, you must specify the x, y, and grouping variables.
#' Optionally, you can pass variables (...) into the labs() in ggplot.
#'
#' @param object model object
#' @param newdata data set to use for prediction (must include all predictors used in model object)
#' @param x_axis x-axis for ggplot
#' @param x_range if choosing "stata" for which_smooth, sequance of x-values for which to produce summary Stata-like margins
#' @param y_axis y-axis for ggplot, predicted value by default
#' @param grouping grouping variable for colour and shape of ggplot points (and optionally LOESS curve)
#' @param which_smooth optional parameter to add LOESS curve
#' * "both" = default; adds a LOESS overall and another LOESS grouped by grouping variable
#' * "single" = adds LOESS curve, ungrouped
#' * "by_grouping" = adds LOESS curve, grouped by the grouping variable
#' * any other input will not add a LOESS curve
#' * "stata" will produce Stata-like margins, using the calculate_stats() function
#' @param point_shape optional, point shape if choosing Stata-like margins. Default value is 16.
#' @param ... optional, to pass on to the labs() in ggplot
#' @seealso [ggplot2::labs()] for more detailed instructions
#' @param legend_title Title for your legend, by default the grouping variable name
#' @param to_return optionally, tell make_margins() what to return
#' * "plot" is the default, and returns a ggplot object
#' * "df" returns the parameter newdata with the addition of a "predicted" column
#' * anything else will return only the "predicted column, which can be easily joined to newdata
#' @param translate TRUE/FALSE (FALSE is default). Should the clock residuals be translated back to unstandardised values?
#' @param sd_ratio optional; standard deviation ratio to translate standardised to unstandardised
#' @param add_null_line optional; adds a dashed red line at y = 0 if TRUE
#'
#' @seealso [base::cbind()] for more detailed instructions
#' @returns either the predicted values (with or without newdata) or a plot, as specified using the to_return argument
#' @export
make_margins <- function(object, newdata, x_axis = "", y_axis = "predicted", grouping = "", x_range = seq(-2, 6, 1), which_smooth = "both", point_shape = 16, bandwidth = 0.5, ..., legend_title = "grouping",
                         to_return = "plot", translate = FALSE, sd_ratio = NULL, add_null_line = FALSE) {

  if (paste0(x_axis, grouping) == "" & to_return == "plot") {
    stop("Please provide x_axis, y_axis, and grouping variable for the plot output")
  }

  if (legend_title == "grouping") {
    legend_title = grouping
  }

  # making margins
  df = newdata %>% cbind(
    predicted = predict(object = object, newdata = newdata)
  )
  if (translate != FALSE & is.null(sd_ratio)) {
    stop("Translate column does not exist in the dataset.")
  } else if (translate != FALSE) {
    df = df %>% mutate(
      predicted2 = predicted
    ) %>% mutate(
      predicted = predicted2*(sd_ratio)
    )
  }
  # plot
  if (which_smooth == "stata" & grouping == "") {
    df_summary = do.call(rbind, lapply(x_range, function(w) calculate_stats(df, x_axis, predicted, w, bandwidth)))
    # error bar width
    eb_width = (max(x_range) - min(x_range))/40
    p = ggplot(
      df_summary,
      aes(x = my_x, y = mean)
      ) +
      geom_point(size = 3, shape = point_shape) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = eb_width) +
      geom_line() +
      labs(
        !!! ensyms(...)
      )
  } else if (which_smooth == "stata" & grouping != "") {
    eb_width = (max(x_range) - min(x_range))/40
    number_groups = length(unique(df[[grouping]]))
    if (number_groups == 1) {
      stop("Grouping variable is unique. Remove grouping parameter or ensure that grouping variable has at least two unique values")
    } else {
      df_summary = do.call(rbind, lapply(x_range, function(w) calculate_stats(df[df[[grouping]] == unique(df[[grouping]])[1],], x_axis, predicted, w, bandwidth))) %>% mutate(
        my_group = unique(df[[grouping]])[1]
      )
      for (i in 2:number_groups) {
        df_summary = rbind(
          df_summary,
          do.call(rbind, lapply(x_range, function(w) calculate_stats(df[df[[grouping]] == unique(df[[grouping]])[i],], x_axis, predicted, w, bandwidth))) %>% mutate(
            my_group = unique(df[[grouping]])[i]
            )
          )
      }
      p = ggplot(
        df_summary,
        aes(x = my_x, y = mean)
      ) +
        geom_point(size = 3, aes(shape = my_group, colour = my_group)) +
        geom_errorbar(aes(ymin = lower, ymax = upper, colour = my_group), width = eb_width) +
        geom_line(aes(colour = my_group)) +
        labs(
          !!! ensyms(...)
        ) +
        theme(plot.caption.position = "panel", plot.caption = element_text(hjust = 0)) +
        guides(
          col = guide_legend(title = legend_title),
          shape = guide_legend(title = legend_title)
        )
    }
  } else {
    p = ggplot(
      df,
      aes_string(x = x_axis, y = "predicted")
    ) +
      geom_point(size=2.5, aes_string(colour = grouping, shape = grouping)) +
      scale_shape_manual(values=c(1,2)) +
      labs(
        !!! ensyms(...)
      ) +
      theme(plot.caption.position = "panel", plot.caption = element_text(hjust = 0)) +
      guides(
        col = guide_legend(title = legend_title),
        shape = guide_legend(title = legend_title)
      )
  }

  if (which_smooth == "both" | which_smooth == "single") {
    p = p + geom_smooth()
  }

  if (which_smooth == "both" | which_smooth == "by_grouping") {
    p = p + geom_smooth(aes_string(colour = grouping))
  }

  if (add_null_line == TRUE) {
    p = p + geom_hline(yintercept = 0, col = "darkred", linetype = "dotted")
  }

  if (to_return == "plot") {
    return(p)
  } else if (to_return == "df") {
    return(df)
  } else {
    return(df$predicted)
  }
}

#' Compare linear regression coefficients
#'
#' @description Allows to empirically compare coefficients for the same predictor in a linear regression model
#' either before and after adjustment for covariates (nested models, not independent)
#' or for different samples (can use the same formula, models assumed independent).
#'
#' @param model1 linear regression model object; reduced model if nested
#' @param model2 linear regression model object; full model if nested
#' @param predictor character string name of coefficient for which to compute difference
#' @param predictor2 if in model2 the coefficient of interest is named something else; NULL by default
#' @param nested TRUE by default
#' @param ci outputs 95\% confidence interval if TRUE
#' @param pval outputs p-value if TRUE
#'
#' @return Returns a vector, which includes difference of coefficients, standard error, and, optionally, 95\%CI, p-value.
#'
#' @seealso
#'     Statistical Methods for Comparing Regression Coefficients Between Models
#'   Author(s): Clifford C. Clogg, Eva Petkova and Adamantios Haritou
#'   Source: American Journal of Sociology, Vol. 100, No. 5 (Mar., 1995), pp. 1261-1293
#'   Published by: The University of Chicago PressStable URL: https://www.jstor.org/stable/2782277
#'
#' @export
coeff_comp = function(model1, model2, predictor, predictor2 = NULL, nested = TRUE, ci = TRUE, pval = TRUE) {

  # model1 is the reduced mode
  # model2 is the full model
  model1t = broom::tidy(model1) %>% filter(term == predictor | term == paste0(predictor, "1"))
  if (is.null(predictor2)) {
    model2t = broom::tidy(model2) %>% filter(term == predictor | term == paste0(predictor, "1"))
  } else {
    model2t = broom::tidy(model2) %>% filter(term == predictor2 | term == paste0(predictor2, "1"))
  }

  # mean-squared errors
  mse1 = mean(model1$residuals^2)
  mse2 = mean(model2$residuals^2)

  d = model1t$estimate - model2t$estimate # difference
  d.sd = ifelse(
    nested == TRUE,
    sqrt(
      model2t$std.error^2 - (model1t$std.error*mse2/mse1)^2
    ), # standard error for difference of nested models
    sqrt(
      model2t$std.error^2 + model1t$std.error^2
    ) # standard error for difference of independent models
  )

  # 95% confidence interval
  d_lower = d - 1.96*d.sd
  d_upper = d + 1.96*d.sd

  my_stat = d/d.sd # t- or z-test, to create p-value
  pvalue = 2*min(
    pnorm(abs(my_stat), lower.tail = TRUE),
    pnorm(abs(my_stat), lower.tail = FALSE)
  )

  if (ci == TRUE & pval == TRUE) {
    output = c("d" = d, "std.error" = d.sd, "ci95.lower" = d_lower, "ci95.upper" = d_upper, "p.val" = pvalue)
  } else if (ci == TRUE & pval == FALSE) {
    output = c("d" = d, "std.error" = d.sd, "ci95.lower" = d_lower, "ci95.upper" = d_upper)
  } else if (pval == TRUE) {
    output = c("d" = d, "std.error" = d.sd, "p.val" = pvalue)
  } else {
    output = c("d" = d, "std.error" = d.sd)
  }

  return(output)
}

#' Compare linear regression models using ANOVA
#'
#' @description Allows to compare up to three linear regression models. Returns the reduced model if ANOVA is statistically significant or the full model. Takes up to three models as input
#'
#' @param model1 linear regression model object; full model
#' @param model2 linear regression model object; reduced model
#' @param model3 optional; linear regression model object; most reduced model
#'
#' @return Returns optimal linear regression model object.
#'
#' @export
return.fm = function(model1, model2, model3 = NULL) {
  # test first round of models (note: model 1 is full, as numbers increase, they get reduced)
  if (paste(as.character(formula(model1)), collapse=" ")==paste(as.character(formula(model2)), collapse=" ")) {
    opt_m = model1
  } else {
    a1 = anova(model1, model2)
    opt_m = if(a1$`Pr(>F)`[2] < 0.05) { model1 } else { model2 }
  }
  if (is.null(model3)) {
    return(opt_m)
  } else if (paste(as.character(formula(opt_m)), collapse=" ")==paste(as.character(formula(model3)), collapse=" ")) {
    return(opt_m)
  } else {
    a2 = anova(opt_m, model3)
    opt_m2 = if(a2$`Pr(>F)`[2] < 0.05) { opt_m } else { model3 }
    return(opt_m2)
  }
}

# this is just for this document; makes three plots for variables that have different distributions by sex
makes3margins <- function(objects, newdata, x_axis, y_axis, x_range, offset = 0, bandwidth, title, x_lab, y_lab, translate = FALSE, sd_ratio = NULL, my_legend_guide = "Sex", my_labels = c("Female", "Male")) { # offset is by how much to offset male from female x_range if applicable; object is list(all, f, m)
  # creating dataframe for non-sex-stratified
  df1 = make_margins(object = objects[[1]], newdata = newdata, x_axis = x_axis, y_axis = y_axis, grouping = "sex", to_return = "df", translate = translate, sd_ratio = sd_ratio) %>%
    dplyr::select(predicted, paste0("z.", x_axis), x_axis, sex)
  # getting stats
  df1.stats <- rbind(
    do.call(rbind, lapply(x_range, function(x) calculate_stats(df1[df1$sex == 0,], x_axis, "predicted", x, bandwidth))) %>% mutate(sex = 0),
    do.call(rbind, lapply(x_range + offset, function(x) calculate_stats(df1[df1$sex == 1,], x_axis, "predicted", x, bandwidth))) %>% mutate(sex = 1)
  ) %>% mutate(sex = factor(sex))

  # plot non-sex-stratified
  eb_width = (max(x_range) - min(x_range))/40 # error bar width
  p.all = ggplot(df1.stats, aes(x = my_x, y = mean)) +
    geom_point(aes(color = sex, shape = sex), size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = sex), width = eb_width) +
    geom_line(aes(color = sex)) +
    labs(
      title = title, x = x_lab, y = y_lab
    ) +
    theme(plot.caption.position = "panel", plot.caption = element_text(hjust = 0)) +
    guides(
      col = guide_legend(title = my_legend_guide),
      shape = guide_legend(title = my_legend_guide)
    ) + scale_shape_discrete(c(15, 16), labels = my_labels) + scale_colour_discrete(labels = my_labels) +
    geom_hline(yintercept = 0, col = "darkred", linetype = "dotted")

  # female plot
  p.f = make_margins(
    object = objects[[2]],
    newdata = newdata[newdata$sex == 0,],
    x_axis = x_axis, y_axis = y_axis, x_range = x_range, bandwidth = bandwidth,
    which_smooth = "stata",
    title = "Females", x = !!enquo(x_lab), y = !!enquo(y_lab), legend_title = "grouping",
    to_return = "plot", add_null_line = TRUE, translate = translate, sd_ratio = sd_ratio
  )

  # male plot
  p.m = make_margins(
    object = objects[[3]],
    newdata = newdata[newdata$sex == 1,],
    x_axis = x_axis, y_axis = y_axis, x_range = x_range + offset, bandwidth = bandwidth,
    which_smooth = "stata",
    title = "Males", x = !!enquo(x_lab), y = !!enquo(y_lab), legend_title = "grouping",
    to_return = "plot", add_null_line = TRUE, translate = translate, sd_ratio = sd_ratio
  )

  return(list(p.all, p.f, p.m))
}

# patchwork function to arrange plots from make3plots
collect_plots = function(list1, list2 = NULL, ylims = TRUE) { # y_lims = TRUE tells it to compute y_lims
  if (ylims == TRUE) {
    if (is.null(list2)) {
      y_lims = lapply(list1, function(x) layer_scales(x)$y$range$range) %>% unlist() %>% range()
      list1 = lapply(list1, function(x) x + ylim(y_lims))
    } else {
      y_lims = lapply(append(list1, list2), function(x) layer_scales(x)$y$range$range) %>% unlist() %>% range()
      list1 = lapply(list1, function(x) x + ylim(y_lims))
      list2 = lapply(list2, function(x) x + ylim(y_lims))
    }
  }
  if (is.null(list2)) {
    return(list1[[1]] / list1[[2]] + list1[[3]] + plot_layout(design = c(area(1,1,8,6), area(9,1,12,3), area(9,4,12,6))))
  } else {
    return((list1[[1]] + list2[[1]] + plot_layout(nrow = 1, guides = "collect")) / (list1[[2]] + list1[[3]] + list2[[2]] + list2[[3]] + plot_layout(nrow = 1)))
  }
}

collect_plots2 = function(p1, p2, p3) {
  return(
    (p1 + guides(col = FALSE, shape = FALSE)) + (p1 %>% get_legend()) + p2 + p3 + plot_layout(
      design = "
  AAB#
  CCDD
  ")
  )
}


#' Make csv into publication-ready excel table
#'
#' @description Make csv into publication-ready excel table
#'
#' @param file.path path of folder in which file of interest is located
#' @param file.name name of csv file (including ".csv")
#'
#' @return Returns a dataframe.
#'
#' @export
csv.to.pub = function(file.path, file.name) {
  file = read.csv(paste0(file.path, "/", file.name)) %>% mutate(
    term = case_when(
      grepl("bmi$", term) ~ "BMI",
      grepl("bmi × ", term) ~ gsub(".*bmi × ", "BMI × ", term),

      grepl("waist$", term) ~ "WC",
      grepl("waist × ", term) ~ gsub(".*waist × ", "WC × ", term),

      grepl("whtr$", term) ~ "WHtR",
      grepl("whtr × ", term) ~ gsub(".*whtr × ", "WHtR × ", term),
      term == "cd8t" ~ "CD8T cells",
      term == "cd4t" ~ "CD4T cells",
      term == "nk" ~ "Natural killers",
      term == "bcell" ~ "B cells",
      term == "mono" ~ "Monocytes",
      term == "gran" ~ "Granulocytes",
      term == "cursmok" ~ "Smoker",
      term == "hai" ~ "Mean lifetime wealth index",
      term == "num.preg.comp" ~ "Number completed pregnancies",
      term == "totmonthsbf" ~ "Total months breastfeeding completed",
      term == "num.preg.comp × totmonthsbf" ~ "Completed pregnancies * Months breastfeeding",
      term == "uw" ~ "Underweight (uw)",
      term == "sex1" ~ "Male sex",
      TRUE ~ term
    )
  )
  return(file)
}


#' Translate the name of the model to something publication-readable
#'
#' @description Make sheet names for the final excel file
#'
#' @param file.name name of csv file (including ".csv")
#'
#' @return Returns a string
#'
#' @export
make.sheet.name = function(file.name) {
  # getting exposure
  p1 = exposures_pub[which(
    grepl(
      gsub("_.*", "", file.name), exposures_pub
    )
  )]

  p2 = ifelse(
    grepl("_pc.csv", file.name),
    paste("PC", mvreg_response_labels[which(
      grepl(
        gsub(".*_", "", gsub("_pc.csv", "", file.name)), mvreg_response_pc
      )
    )]),
    mvreg_response_labels[which(
      grepl(
        gsub(".*_", "", gsub(".csv", "", file.name)), mvreg_response_nonpc
      )
    )]
  )

  sheetname = paste(p1, p2)

  return(sheetname)

}


#' Put dataframes together into an excel file
#'
#' @description put a list of .csv files into a .xlsx files, where each tab (sheet) is named after the csv
#'
#' @param file.list list of .csv files (no path, just file names)
#' @param file.path path to the folder containing the list of .csv files to be combined
#' @param path.to.write.to path of the folder to which you want to write out the final .xlsx file
#' @param file.name name of .xlsx file to be made and writen out
#'
#' @return Returns a .xlsx file
#'
#' @export
tabble = function(file.list, file.path, path.to.write.to, file.name) {
  n = length(file.list) # for loop
  excel.name = paste0(path.to.write.to, "/", file.name, ".xlsx") # name of new Excel file (full path)

  # first file to write
  # write.xlsx(
  #   x = csv.to.pub(file.path, file.list[1]),
  #   file = excel.name,
  #   sheetName = sheet.name
  #   )

  workbook = createWorkbook()

  for (i in 1:n) {
    addWorksheet(
      workbook,
      gsub(".csv", "", my.files[i])
    )
    writeData(
      wb = workbook,
      sheet = gsub(".csv", "", my.files[i]),
      x = csv.to.pub(file.path, file.list[i])
    )
  }

  saveWorkbook(
    wb = workbook,
    file = excel.name
  )

  return(workbook)

}

#' Print out the coefficient comparison table
#'
#' @description print out the coefficient comparison table
#'
#' @return Returns a kable
#'
#' @export
# po = print out
po.kable = function(tibble, caption = " ") {
  return(
    tibble %>% as.data.frame() %>% mutate(across(1:5, as.numeric)) %>% rename(
      "Difference" = 1,
      "Std.Error" = 2,
      "p.value" = 5,
      "Compared" = 6
    ) %>% mutate(
      `95% CI` = paste0("(", sprintf("%.3f", ci95.lower), ", ", sprintf("%.3f", ci95.upper), ")")
    ) %>% dplyr::select(
      "Compared", "Difference", "95% CI", "Std.Error", "p.value"
    ) %>% kable(caption = caption, digits = 3) %>% kable_styling("striped", full_width = F)
  )
}
