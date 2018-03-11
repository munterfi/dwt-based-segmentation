#' Extract movement parameters from track
#'
#' @param xCoord longitude of track points in decimal degrees
#' @param yCoord Latitude of track points in decimal degrees
#' @param tCoord time lag between track points (not cumsum!)
#'
#' @return Matrix of movement parameters
#' @export
#'
#' @examples
#' movementParameters(xCoord, yCoord, tCoord)
movementParameters <- function(xCoord, yCoord, tCoord) {
  xyCoord <- as.matrix(cbind(xCoord, yCoord))
  # empty arrays
  sinuosityValues = numeric(nrow(xyCoord))
  tortuosityValues = numeric(nrow(xyCoord))
  velocityValues = numeric(nrow(xyCoord))
  accelerationValues = numeric(nrow(xyCoord))
  turnAngleValues = numeric(nrow(xyCoord))
  angularVelocityValues = numeric(nrow(xyCoord))
  meanderValues = numeric(nrow(xyCoord))
  # loop through track
  w=1
  for (i in (1+w):(nrow(xyCoord)-1-w)) {
    timeStep=tCoord[i]
    points = rbind(xyCoord[i-w, ], xyCoord[i, ], xyCoord[i+w, ])
    deltaX12 = xyCoord[i,1]-xyCoord[i-w,1]
    deltaY12 = xyCoord[i,2]-xyCoord[i-w,2]
    deltaX23 = xyCoord[i+w,1]-xyCoord[i,1]
    deltaY23 = xyCoord[i+w,2]-xyCoord[i,2]
    dist = c(dist(points))
    # sinuosity
    distanceTraveled = 0
    for (q in (i-w):(i+w-1)) {
      pointsW = rbind(xyCoord[q, ], xyCoord[q+1, ])
      distanceW = c(dist(pointsW))
      distanceTraveled = distanceTraveled + distanceW
    }
    endpoints = rbind(xyCoord[i-w, ], xyCoord[i+w, ])
    actualDistance = c(dist(endpoints))
    sinuosityValues[i] = distanceTraveled/actualDistance
    # velocity
    v1 = sqrt((deltaX12)^2+(deltaY12)^2)/(timeStep*w)
    v2 = sqrt((deltaX23)^2+(deltaY23)^2)/(timeStep*w)
    velocityValues[i]=(v1+v2)/2
    # acceleration
    accelerationValues[i] = (velocityValues[i]-velocityValues[i-w])/(timeStep*w)
    # turning angle
    turnAngleValues[i] = -(atan2(deltaX12*deltaY23-deltaX23*deltaY12,deltaX12*deltaX23+deltaY12*deltaY23))
    # meander values
    if (dist[2] != 0){
      meanderValues[i] = turnAngleValues[i]/dist[2]
    } else {
      meanderValues[i] = 0
    }
    # angular velocity values
    angularVelocityValues[i] = turnAngleValues[i]/(timeStep*w);
    # tortuosity values
    tortuosityValues[i] = abs(velocityValues[i])*cos(turnAngleValues[i]);
  }
  MPvel = velocityValues
  MPsin = sinuosityValues
  MPtor = tortuosityValues
  MPacc = accelerationValues
  MPturn = turnAngleValues
  MPang = angularVelocityValues
  MPmeander = meanderValues
  return(list(MPvel = velocityValues, MPsin = sinuosityValues, MPtor = tortuosityValues,
         MPacc = accelerationValues, MPturn = turnAngleValues, MPang = angularVelocityValues,
         MPmeander = meanderValues))
}
