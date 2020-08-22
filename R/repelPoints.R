# Repel points (shift points 40% to their inter cluster distance.)
shiftPoints <- function(clust, leftDist, rightDist) {
  if (nrow(clust) == 1) {
    clust$pos2 <- clust$pos
    return(clust)
  } else if (nrow(clust) == 2) {
    # if distance is huge between right and left clusters, bring it down to 50.
    leftDist <- ifelse(test = leftDist > 100, yes = 50, no = leftDist)
    rightDist <- ifelse(test = rightDist > 100, yes = 50, no = rightDist)

    startPoint <- clust[1, pos] - as.integer(leftDist * 0.2)
    endPoint <- clust[nrow(clust), pos] + as.integer(rightDist * 0.2)
    clust$pos2 <- c(startPoint, endPoint)
    return(clust)
  } else {
    leftDist <- ifelse(test = leftDist > 100, yes = 50, no = leftDist)
    rightDist <- ifelse(test = rightDist > 100, yes = 50, no = rightDist)

    startPoint <- clust[1, pos] - as.integer(leftDist * 0.4)
    endPoint <- clust[nrow(clust), pos] + as.integer(rightDist * 0.4)
    posPoints <- clust[2:(nrow(clust) - 1), pos]

    dist <- rbind(abs(posPoints - startPoint), abs(posPoints - endPoint))
    dist <- sapply(apply(dist, 2, function(x) which(x == max(x))), "[[", 1)


    tempStart <- startPoint
    tempEnd <- endPoint

    for (i in 1:length(posPoints)) {
      if (dist[i] == 1) {
        posPoints[i] <- tempStart + as.integer(0.4 * abs(tempStart - posPoints[i]))
        tempStart <- posPoints[i]
      } else {
        posPoints[i] <- tempEnd - as.integer(0.4 * abs(tempEnd - posPoints[i]))
        tempEnd <- posPoints[i]
      }
    }
    clust$pos2 <- c(startPoint, posPoints, endPoint)
    return(clust)
  }
}

# cluster
repelPoints <- function(dat, protLen, clustSize) {
  dat <- data.table(dat)

  dat <- dat[order(dat$pos), ]
  dat$distance <- c(0, diff(dat$pos))

  ## cluster identification
  mergeDist <- clustSize # hard coded inter event distance.
  cdf <- c()
  cluster <- 1
  n <- 1
  start <- end <- dat[1, pos]

  for (i in 2:nrow(dat)) {
    dist <- dat[i, distance]
    position <- dat[i, pos]

    if (dist <= mergeDist) {
      end <- position
      n <- n + 1
    } else {
      tempdf <- data.frame(cluster = paste("cluster", cluster, sep = "_"), start = start, end = end, n = n)
      n <- 1
      cdf <- rbind(cdf, tempdf)
      cluster <- cluster + 1
      start <- end <- position
    }
  }
  cdf <- rbind(cdf, data.frame(cluster = paste("cluster", cluster, sep = "_"), start = start, end = end, n = n))

  # Assign clusters to main table
  cdf2 <- c()
  for (i in 1:nrow(cdf)) {
    start <- cdf[i, "start"]
    end <- cdf[i, "end"]
    tdat <- dat[pos >= start & pos <= end]
    tdat[, cluster := paste("cluster", i, sep = "_")]
    cdf2 <- rbind(cdf2, tdat)
  }

  # Cluster summary and inter cluster distance
  clustSum <- cdf2[, .(min = min(pos), max = max(pos)), by = cluster]
  clustEnd <- clustSum[1, max]
  clustDist <- 0
  if (nrow(clustSum) > 1) {
    for (i in 2:nrow(clustSum)) {
      clustDist <- c(clustDist, clustSum[i, min] - clustEnd)
      clustEnd <- clustSum[i, max]
    }
  }

  clustSum[, dist := clustDist]

  ## repel points
  cdf3 <- c()
  clusters <- clustSum[, cluster]

  for (i in 1:length(clusters)) {
    tempClust <- cdf2[cluster == clusters[i]]
    # print(tempClust)
    if (i == 1) {
      cdf3 <- rbind(cdf3, shiftPoints(clust = tempClust, leftDist = tempClust[1, pos], rightDist = clustSum[2, dist]))
      # In case repelled points cross protein length, bring it down to protein len.
      cdf3$pos2 <- ifelse(test = cdf3$pos2 >= protLen, yes = protLen, no = cdf3$pos2)
    } else if (i == length(clusters)) {
      cdf3 <- rbind(cdf3, shiftPoints(clust = tempClust, leftDist = clustSum[i, dist], rightDist = protLen))
      cdf3$pos2 <- ifelse(test = cdf3$pos2 >= protLen, yes = protLen, no = cdf3$pos2)
    } else {
      cdf3 <- rbind(cdf3, shiftPoints(clust = tempClust, leftDist = clustSum[i, dist], rightDist = clustSum[i + 1, dist]))
      cdf3$pos2 <- ifelse(test = cdf3$pos2 >= protLen, yes = protLen, no = cdf3$pos2)
    }
  }

  return(cdf3)
}
