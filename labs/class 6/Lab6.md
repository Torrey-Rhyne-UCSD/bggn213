Lab 6: R Functions
================
Torrey (A14397504)

## All about functions in R

Every function in R has at least 3 things: - name (you choose) -
arguments (the input(s)) - the body

Today we will write a function to grade a class of student assignment
scores.

Q1. Write a function grade() to determine an overall grade from a vector
of student homework assignment scores dropping the lowest single score.
If a student misses a homework (i.e. has an NA value) this can be used
as a score to be potentially dropped. Your final function should be
adequately explained with code comments and be able to work on an
example class gradebook such as this one in CSV format:
“https://tinyurl.com/gradeinput” \[3pts\]

``` r
# load gradebook
gradebook <- read.csv("https://tinyurl.com/gradeinput", row.names = 1)
# load practice data
student1 <- c(100, 100, 100, 100, 100, 100, 100, 90)
student2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
student3 <- c(90, NA, NA, NA, NA, NA, NA, NA)

# write function
grade <- function(x) {
  # convert NAs to 0
  x[is.na(x)] <- 0 
  # drop the lowest score
  x_dropped <- x[-which.min(x)]
  # calculate the average
  mean(x_dropped)
}

# grade practice students
grade(student1)
```

    [1] 100

``` r
grade(student2)
```

    [1] 91

``` r
grade(student3)
```

    [1] 12.85714

Q2. Using your grade() function and the supplied gradebook, Who is the
top scoring student overall in the gradebook? \[3pts\]

``` r
# who got the highest score?
which.max(apply(gradebook,1,grade))
```

    student-18 
            18 

``` r
# what is their score?
max(apply(gradebook,1,grade))
```

    [1] 94.5

Q3. From your analysis of the gradebook, which homework was toughest on
students (i.e. obtained the lowest scores overall? \[2pts\]

``` r
# Answer depends on interpretation. 
# Not doing the HW (NA) probably means the student was busy, not that it was tough. 
# I will choose to remove NAs.

which.min(apply(gradebook, 2, mean, na.rm = T))
```

    hw3 
      3 

``` r
min(apply(gradebook, 2, mean, na.rm = T))
```

    [1] 80.8

Q4. Optional Extension: From your analysis of the gradebook, which
homework was most predictive of overall score (i.e. highest correlation
with average grade score)? \[1pt\]

``` r
# store final score for each student
scores <- as.data.frame(apply(gradebook,1,grade)); colnames(scores) = "scores"
# add to gradebook 
gradebook <- cbind(gradebook,scores)

# calculate correlation between each HW and final score
# again, choosing to remove NA values
cor <- cor(gradebook[,1:5], gradebook[,6], method = "pearson", use = "pairwise.complete.obs")

# which HW is most predictive of the final score? what is the correlation?
which.max(cor); max(cor)
```

    [1] 2

    [1] 0.6114277

Q5. Make sure you save your Quarto document and can click the “Render”
(or Rmarkdown”Knit”) button to generate a PDF foramt report without
errors. Finally, submit your PDF to gradescope. \[1pt\]
