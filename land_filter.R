# Build and test case code for the land filter dynamic library
require(BuildSys)

Project <- new("BSysProject", "./", CXXFLAGS="-D__R_MODULE__", Debug=T)
make(Project, "clean")
make(Project)
loadLibrary(Project)

vcDebug(Project)

Context <- .External('Initialise', "./coastline/GSHHS_f_L1.shp")

require(ozmaps)
xlim <- c(145,148)
ylim <- c(-43.7,-42)
ozmap(xlim=xlim, ylim=ylim)

max_speed       <- as.double(8.0)
seconds_in_days <- 3600 * 24 
dt              <- 1

start <- c(145,-43.0)
finish <- c(146,-43.0)

res <- .External('SpeedAndLandLimit', Context, start[1], start[2], finish[1], finish[2], dt * seconds_in_days, max_speed)
lines(c(start[1],res[1]), c(start[2],res[2]), col="blue")

start <- c(146,-43.0)
finish <- c(145,-43.0)

res <- .External('SpeedAndLandLimit', Context, start[1], start[2], finish[1], finish[2], dt * seconds_in_days, max_speed)
lines(c(start[1],res[1]), c(start[2],res[2]), col="green")

start <- c(149,-43.0)
finish <- c(145,-43.0)

res <- .External('SpeedAndLandLimit', Context, start[1], start[2], finish[1], finish[2], dt * seconds_in_days, max_speed)
lines(c(start[1],res[1]), c(start[2],res[2]), col="purple")

start <- c(149,-43.0)
finish <- c(146,-43.0)

res <- .External('SpeedAndLandLimit', Context, start[1], start[2], finish[1], finish[2], dt * seconds_in_days, max_speed)
lines(c(start[1],res[1]), c(start[2],res[2]), col="purple")

start <- c(147,-43.7)
finish <- c(147,-42.0)

res <- .External('SpeedAndLandLimit', Context, start[1], start[2], finish[1], finish[2], dt * seconds_in_days, max_speed)
lines(c(start[1],res[1]), c(start[2],res[2]), col="orange")

start <- c(147,-42.0)
finish <- c(147,-43.7)

res <- .External('SpeedAndLandLimit', Context, start[1], start[2], finish[1], finish[2], dt * seconds_in_days, max_speed)
lines(c(start[1],res[1]), c(start[2],res[2]), col="grey")
