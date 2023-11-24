# This is just a text file.

# base R graphics
x <- 1:50
plot(x) # value vs. index
plot(x, sin(x)) # now x and y, but this is ugly - argue
plot(x, sin(x), type = "l")
plot(x, sin(x), type = "b", col = "pink", lwd = 5)

# use packages instead of base R
# for plots - ggplot


