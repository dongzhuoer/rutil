# library(ggplot2)



# ## to do: use list instead of eval(parse)


# plot_blast <- function(json.file, fields = c( 'bit_score', 'score', 'evalue', 'identity', 'align_len', 'gaps')) {
#     json <- jsonlite::read_json(json.file);

#     select <- function(field) {
#         sapply(json$BlastOutput2$report$results$search$hits, function(x){x$hsps[[1]][[field]]})
#     }

#     eval(parse(text = paste0('df <- data.frame(x=1:100', paste0(', ', fields, '=select("', fields, '")', collapse = ''), ')')))

#     g <- eval(parse(text = paste0('cowplot::plot_grid(', paste0('ggplot2::ggplot(df)+ggplot2::geom_line(aes(x=x,y=', fields, '))', collapse = ','), ')')))
#     print(g)
#     g
# }

# lapply(dir('../blast', pattern = 'json', full.names = T), plot_blast)




