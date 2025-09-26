shh <- suppressPackageStartupMessages
shh(library(tidyverse))
shh(library(data.table))

options(repr.matrix.max.cols=50, repr.matrix.max.rows=100)
options(dplyr.summarise.inform = FALSE)
options(warn = -1)

mn <- function(x) round(mean(x, na.rm = TRUE), 2)
corner <- function(mat, k = 5) mat[1:k,1:k]

### dplyr ###
gb <- group_by
ug <- ungroup
ar <- arrange
mu <- mutate
ij <- inner_join
lj <- left_join
rj <- right_join
se <- select
tm <- transmute
su <- summarise
rw <- rowwise
re <- reframe
fi <- filter
pu <- pull
sp <- spread
ma <- mutate_all
ren <- rename_with
uni <- unique
ga <- gather
m <- grepl

ds <- function(w = 10, h = 6){
  options(repr.plot.width = w, repr.plot.height = h)
} 
df <- data.frame
ac <- function(i) as.character(i)

top <- function(data, group_col) {
  data %>%
   group_by({{ group_col }}) %>%
   summarise(ct := n(), .groups = "drop") %>%
   arrange(desc(ct))
}
