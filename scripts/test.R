library(dplyr)
library(magrittr)
t = tibble(a = c(1,3), b = c(2,4))
t


t = t %>% mutate(ten = 10)
t

my_var = "columnname"
t = t %>% mutate({{my_var}} := ten + 1)
t
my_var2 = "columnname2"
t %>% mutate({{my_var2}} := .data[[my_var]] +1)
t %>% mutate({{my_var2}} := .data[[my_var]] +1)