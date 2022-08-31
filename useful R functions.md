# Generic R functions

## Priceless advice

[What They Forgot to Teach You About R](https://whattheyforgot.org/) 

- Nominare i file evitando spazi, accenti e puntuazione. Usa “-” per la human redability e “_” per isolare i metadata dal nome (es. "2020-02-10_sample1_WB") e identificarli facilmente con regex: 

```r
metadata <- stringr::str_split_fixed(x, [_\\.]
```

* Crea Rprojects: avviando un notebook da progetto sei sicuro che la wd sia settata nella cartella di appartenenza e che l’ambiente delle variabili sia quello giusto. Non usare mai patwhay assoluti, soprattutto mai usare setwd(). Lascia che sia l’Rproj a richiamare la WD

* Riavviare R prima di ogni esecuzione perchè:
  
  * potresti star usando funzioni omonime di librerie precedentemente caricate
  
  * un altro progetto potrebbe averti portato nella WD sbagliata
- lasciare settato options(stringsAsFactors = FALSE)

- Carica la libreria tidyverse per ultima, in modo da evitare che altre librerie mascherino le sue funzioni

## Basic R

#### define a function

```r
myfunction <- function(arg1, arg2, ... ){                  
  # statements  
  return(object)  
}
```

#### match operator

```r
x %in% c("value1", "value2", "value3")
```

#### opposite of match operator

```r
'%!in%' <- function(x,y)!('%in%'(x,y))
```

## Magrittr

#### Pipe operator

```r
f(x) # is equivalent to
x %>% f

f(x, parameters) # is equivalent to
x %>% f(parameters) 

h(g(f(x))) # is equivalent to
x %>% f %>% g %>% h 

f(y, x) # is equivalent to
x %>% f(y, .)

f(y, data = x) # is equivalent to
x %>% f(y, data = .)

f(g(x))
x %>% {f(g(.))}
```

#### Compound assigmnent

```r
data$some_variable <- data$some_variable %>% transform 
# is equivalent to
data$some_variable %<>% transform
```

## Data exploration

```r
summary(iris)
skimr::skim(iris)
```

## Dplyr

#### Matches at least with one element of the list

```r
3 %in% c(1,2,3,4,5,6,7,8)  
[1] TRUE
```

#### Select rows by criteria

```r
flights[flights$month == 1 & flights$day == 1, ]
# is equivalent to
flights %>% filter(month == 1, day == 1)
# multiple arguments are equivalent to &
flights %>% filter(month == 1 & day == 1)
```

#### Select only certains columns

```r
# select a single column

flights %>% select(tailnum)

# select multiple column

flights %>% select(tailnum, day, month)

# select and rename the column

flights %>% select(tail_num = tailnum)

# select all but some columns

flights %>% select(-(year:day)) # or
flights %>% select(-c(year, month))

# select column based on theire name
flights %>% select(starts_with(c("one", "th")))
```

#### Rename a column

```r
flights %>% rename(tail_num = tailnum)
```

#### Move a colon to the first position

```r
df <- df %>% select(fisrt_column, second_column, everything())
```

#### Create a new column with formula

```r
mutate(flights,  
 gain = arr_delay - dep_delay,  
 gain_per_hour = gain / (air_time / 60)  
)
# use the function transmutate() if you want to remove the used columns
```

#### Replace a single value by criterion

```r
df <- df %>% mutate(height = replace(height, name == “Mike”, NA))
```

#### Summarize

```r
summarise(flights,  
  delay.mean = mean(dep_delay, na.rm = TRUE), delay.sd = sd(dep_delay, na.rm = TRUE)  
)
```

## tidyr

#### gather() - wide to long format

![](https://lh5.googleusercontent.com/sW47p8p8ifQo4ArWuXIWriSrw8d_7VPL94oa_OW2bbllXmjsxsB6yFE57Z9Y91isS2OoEbJiNvld_Hdm0qx0vhoftAqs_dGRawPhQV0L7UX0V-NRr1rIw5znSeLx_h92TseLCjsU)

#### spread() - long to wide format![](https://lh3.googleusercontent.com/paaS7OyTGxSAE5z9uOil23kgut2mBSF7qFCy8Uqx0slVcqnI2QQ-fsRAvGiQjtH8K6xFhJ75lF2JkgQOsO-dJmvFgBRQj39o79paDOvsEwJ2sJcv7uVXvmGs9BjcFZBOesFfnXMs)

### transpose a dataframe 
```r
df <- df %>%
  rownames_to_column() %>%
  gather(var, value, -rowname) %>% 
  spread(rowname, value) 
```

## Favorite structures

#### Rename a specific column

```r
# rename the second column
colnames(data.frame)[2] <- "newname2"
```

#### Select all the rows that contain a certain match inside a specified column

```r
# selects all the rows containing a certain match inside a certain col
data.frame %>% filter(str_detect(col, fixed("match")))
```

#### Split the string inside a column into multiple column

```r
# will split a string like "CTRL-2day LowGlucose"
data.frame %>% separate(col, # the column to split  
 c("col1", "col2", "col3"), # names of the columns  
 sep = "[:punct:]|[:blank:]", # separator  
 remove = TRUE, # remove the original column  
)
```

#### Change Factor levels

```r
data.frame$col %<>% factor(levels = c("first", "second", "third"))
```

#### Calculate Standard Error

```r
st.err <- function(x) {sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))}
```

## ggplot2

#### Convert Y axis in percentage

```r
scale_y_continuous(labels = function(x) paste0(x*100, "%"))
```

#### Rename the legend

```r
+ scale_fill_x(name = 'Title')
```

#### Remove legend for one aesthetics

```r
+  guides(fill = FALSE)
# or
+ scale_fill_x(guide = FALSE)
```

#### Declare the specific color for each factor level

```r
myColors <- c("level1" = "#000000", "level2" = "#003300") 
# or
library(RColorBrewer)  
myColors  <- brewer.pal(length(levels(df$factors)), "Set1") 

names(myColors) <- levels(df$factors) 
+ scale_colour_manual(values = myColors)
```

## Generate random data

```r
set.seed(999)  
n = 1000  
df = data.frame(factors = sample(letters[1:8], n, replace = TRUE),  
x = rnorm(n), y = runif(n))

x <- runif(N)
y <- 5 * x + 3 + rnorm(N)

```
