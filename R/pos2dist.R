#' Distance between two points on earth's surface
#'
#' Calculates distance between two points on earth's surface
#' given by their latitude-longitude pair.
#' Input lag1,lon1,lag2,lon2 are in degrees, without 'NSWE' indicators.
#' Input method is 1 or 2. Default is 1.
#' Method 1 uses plane approximation,
#' only for points within several tens of kilometers (angles in rads):
#' \code{d = sqrt(R_equator^2*(lag1-lag2)^2 + R_polar^2*(lon1-lon2)^2*cos((lag1+lag2)/2)^2)}
#' Method 2 calculates sphereic geodesic distance for points farther apart,
#' but ignores flattening of the earth:
#' \code{d = R_aver * acos(cos(lag1)cos(lag2)cos(lon1-lon2)+sin(lag1)sin(lag2))}
#' Flora Sun, University of Toronto, Jun 12, 2004.
#'
#' @param lag1 Numeric skalar in degrees
#' @param lon1 Numeric skalar in degrees
#' @param lag2 Numeric skalar in degrees
#' @param lon2 Numeric skalar in degrees
#' @param method by default set to 1 (plane approximation), other values use spheric geodesic distance
#'
#' @return Dist in km, or -99999 if input argument(s) is/are incorrect.
#' @export
#' 
#' @note The original MatLab version on this function stems from the Soleymani et al. (2017).
#' @author M. Unterfinger (Department of Geography, University of Zurich)
#'
#' @examples
#' pos2dist(lag1,lon1,lag2,lon2,method=1)
pos2dist <- function(lag1,lon1,lag2,lon2,method=1){
  if(missing(lag1,lon1,lag2,lon2)){
    warning('Number of input arguments error! distance = -99999')
    return(dist = -99999)
  }
  if (abs(lag1)>90 | abs(lag2)>90 | abs(lon1)>360 | abs(lon2)>360) {
    warning('Degree(s) illegal! distance = -99999')
    return(dist = -99999)
  }
  if (lon1 < 0) {
    lon1 = lon1 + 360
  }
  if (lon2 < 0) {
    lon2 = lon2 + 360
  }
  # Default method is 1.
  if (method == 1){
    km_per_deg_la = 111.3237
    km_per_deg_lo = 111.1350
    km_la = km_per_deg_la * (lag1-lag2)
    # Always calculate the shorter arc.
    if (abs(lon1-lon2) > 180) {
      dif_lo = abs(lon1-lon2)-180
    } else {
      dif_lo = abs(lon1-lon2)
    }
    km_lo = km_per_deg_lo * dif_lo * cos((lag1+lag2)*pi/360)
    dist = sqrt(km_la^2 + km_lo^2)
  } else {
    R_aver = 6374
    deg2rad = pi/180
    lag1 = lag1 * deg2rad
    lon1 = lon1 * deg2rad
    lag2 = lag2 * deg2rad
    lon2 = lon2 * deg2rad
    dist = R_aver * acos(cos(lag1)*cos(lag2)*cos(lon1-lon2) + sin(lag1)*sin(lag2))
  }
}