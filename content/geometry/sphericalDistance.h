/**
 * Author: Ulf Lundstrom
 * Date: 2009-04-07
 * License: CC0
 * Source: My geometric reasoning
 * Description: Returns the shortest distance on the sphere with radius radius between the points with azimuthal angles (longitude) f1 ($\phi_1$) and f2 ($\phi_2$) from x axis and zenith angles (latitude) t1 ($\theta_1$) and t2 ($\theta_2$) from z axis (0 = north pole). All angles measured in radians. The algorithm starts by converting the spherical coordinates to cartesian coordinates so if that is what you have you can use only the two last rows. d and s are the lengths of the difference and the sum of the two unit vectors; the angle between them is 2 atan2(d, s), which stays accurate (and never NaN) for nearby and antipodal points alike.
 * Status: tested on kattis:airlinehub
 */
#pragma once

double sphericalDistance(double f1, double t1,
		double f2, double t2, double radius) {
	double x1 = sin(t1)*cos(f1), y1 = sin(t1)*sin(f1);
	double x2 = sin(t2)*cos(f2), y2 = sin(t2)*sin(f2);
	double z1 = cos(t1), z2 = cos(t2);
	double dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
	double sx = x2 + x1, sy = y2 + y1, sz = z2 + z1;
	double d = sqrt(dx*dx + dy*dy + dz*dz);
	double s = sqrt(sx*sx + sy*sy + sz*sz);
	return radius*2*atan2(d, s);
}
