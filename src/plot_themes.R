plot_theme = function (font_size = 20) {
    theme_bw(base_size = 12, base_family = 'Helvetica') +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5),
              text = element_text(size = font_size, colour = 'black'),
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.ticks.x = element_line(colour = 'black'),
              axis.ticks.y = element_line(colour = 'black'),
              axis.text.x = element_text(colour = 'black'),
              axis.text.y = element_text(colour = 'black'))
}

plot_theme_facet = function (font_size = 20) {
    theme_bw(base_size = 12, base_family = 'Helvetica') +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5),
              text = element_text(size = font_size, colour = 'black'),
              strip.background = element_blank(),
              panel.border = element_rect(colour="black"),
              axis.ticks.x = element_line(colour = 'black'),
              axis.ticks.y = element_line(colour = 'black'),
              axis.text.x = element_text(colour = 'black'),
              axis.text.y = element_text(colour = 'black'))
}