##==============================================================================
## colorblindPalette.R
##
## Source:
## http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
## http://www.somersault1824.com/wp-content/uploads/2015/02/color-blindness-palette.png
##
## This script is used to generate color combinations that resemble typical 
## colors but are friendlier to those with some form of colorblindness. We also
## recommend downloading the Color Oracle application (https://colororacle.org/)
##
## Code by Tony Wong (aewsma@rit.edu)
##==============================================================================
## Copyright 2019 Tony Wong
## This file is part of GEOCARB-calibration.
## GEOCARB-calibration is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## GEOCARB-calibration is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GEOCARB-calibration.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

mycol=rbind(
              c(0,0,0),
              c(0,73,73),
              c(0,146,146),
              c(255,109,182),
              c(255,182,119),
              c(73,0,146),
              c(0,109,219),
              c(182,109,255),
              c(109,182,255),
              c(182,219,255),
              c(146,0,0),
              c(146,73,0),
              c(219,109,0),
              c(36,255,36),
              c(255,255,109)
            )
mycol=mycol/max(mycol)


##==============================================================================
## End
##==============================================================================
