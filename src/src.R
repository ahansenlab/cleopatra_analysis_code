---
title: "Untitled"
output: html_document
date: "2024-10-24"
---

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

celltypes = c('GM12878', 'HCT116', 'K562', 'H1')
celltype_colours = c('#23344d', '#136f63', '#C64191', '#dd8c22')

loop_order = c('E-P', 'P-P', 'E-E', 'E/P-other', 'CTCF-CTCF', 'CTCF-other', 'other-other')
loop_colours = c('dodgerblue3', 'dodgerblue4', 'skyblue4', '#8C9CB1', '#CC3311', '#a3280d', 'gray')

capture_regions = read_tsv('../../Loci selection/v1_v2_selected_loci.txt',
         col_names = c('locus_id', 'chrom','start', 'end'),
         show_col_types = FALSE) |> 
    mutate(locus_id = locus_id - 1) |> 
    mutate(locus_id = paste0('region', locus_id)) |> 
    makeGRangesFromDataFrame()

capture_regions_200bp_tiles = reduce(tile(capture_regions, width = 200), c)
capture_regions_1kb_tiles = reduce(tile(capture_regions, width = 1000), c)

make_ginteractions_from_df = function(df, extra_cols = NULL){
    
    gr1 = df |> 
        dplyr::select(seqnames1, start1, end1, all_of(extra_cols)) |> 
        makeGRangesFromDataFrame(seqnames = 'seqnames1', start.field = 'start1', end.field = 'end1',
                                 keep.extra.columns = TRUE)
    gr2 = makeGRangesFromDataFrame(df, seqnames = 'seqnames2', start.field = 'start2', end.field = 'end2')
    
    GInteractions(gr1, gr2)
}

read_mcool_files = function(path, re, strip){
    
    fs::dir_ls(path = path, regexp = re) |> 
    as.character() |> 
    purrr::set_names(function(x) {
    str_sub(fs::path_file(x), end = - str_length(strip) - 1)
    }) |> 
    map(~read_tsv(.x, show_col_types = FALSE)) |> 
    map(~select(.x, start1, start2, balanced)) |> 
    map(~filter(.x, !is.na(balanced)))
}

get_anchors_from_interactions = function(df, extend=FALSE){
    
    anchor1_df = df |> 
        dplyr::select(-seqnames2, -start2, -end2) |> 
        dplyr::rename(seqnames = seqnames1, start = start1, end = end1)
    
    anchor2_df = df |> 
        dplyr::select(-seqnames1, -start1, -end1) |> 
        dplyr::rename(seqnames = seqnames2, start = start2, end = end2)
    
    if (extend==FALSE){
        anchors = bind_rows(anchor1_df, anchor2_df) |> 
        distinct(seqnames, start, end) |>
        makeGRangesFromDataFrame()
    } else {
        anchors = bind_rows(anchor1_df, anchor2_df) |> 
        distinct(seqnames, start, end) |>
        mutate(start = start - extend, end = end + extend) |> 
        makeGRangesFromDataFrame()
    }
}

filter_loops_by_loop_class = function(annotated_loops, loop_classes){
    
    annotated_loops |> 
        filter(loop_class %in% loop_classes)
}

write_seqs = function(things, filename, gi=FALSE){
    if (gi == TRUE){
        seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, get_anchors_from_interactions(things))
    } else {
        seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, things)
    }
    names(seqs) = paste0('CRE_anchor_', c(1:length(seqs)))
    writeXStringSet(seqs, filename)
}

read_chip = function(filename){
    read_tsv(filename, col_names = c('chrom', 'start', 'end', paste0('crap', c(1:3)), 'score', 'crap4', 'crap5', 'summit'), show_col_types = FALSE) |> 
    mutate(summit_coord = start + summit) |> 
    select('chrom', 'start', 'end', 'score', 'summit_coord') |> 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) 
}

