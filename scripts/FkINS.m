
This one is in terms of A_xb, A_yb, A_zb etc. 

Fk =
 
[                                                                                                           tan(theta)*(omega_y*cos(phi)-omega_z*sin(phi)),                                                                                                     (1+tan(theta)^2)*(omega_y*sin(phi)+omega_z*cos(phi)),                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0]
[                                                                                                                       -omega_y*sin(phi)-omega_z*cos(phi),                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0]
[                                                                                                           (omega_y*cos(phi)-omega_z*sin(phi))/cos(theta),                                                                                              (omega_y*sin(phi)+omega_z*cos(phi))/cos(theta)^2*sin(theta),                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0]
[                                  (sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi))*conj(A_yb)+(cos(phi)*sin(psi)-sin(phi)*sin(theta)*cos(psi))*conj(A_zb),                                          -sin(theta)*cos(psi)*conj(A_xb)+sin(phi)*cos(theta)*cos(psi)*conj(A_yb)+cos(phi)*cos(theta)*cos(psi)*conj(A_zb), -cos(theta)*sin(psi)*conj(A_xb)+(-cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi))*conj(A_yb)+(sin(phi)*cos(psi)-cos(phi)*sin(theta)*sin(psi))*conj(A_zb),                                                                                                                                                        0,                                                                                                                   (4+lon_dot*sin(lat))*diff(conj(ve),ve),                                                                                                                               -lat_dot*diff(conj(vd),vd),                                                                                                                                lon_dot*cos(lat)*conj(ve),                                                                                                                                                        0,                                                                                                                                                        0]
[                                (-sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi))*conj(A_yb)+(-cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi))*conj(A_zb),                                          -sin(theta)*sin(psi)*conj(A_xb)+sin(phi)*cos(theta)*sin(psi)*conj(A_yb)+cos(phi)*cos(theta)*sin(psi)*conj(A_zb),  cos(theta)*cos(psi)*conj(A_xb)+(-cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi))*conj(A_yb)+(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi))*conj(A_zb),                                                                              (-5380631166442341/36893488147419103232-lon_dot*sin(lat))*diff(conj(vn),vn),                                                                                                                                                        0,                                                                              (-5380631166442341/36893488147419103232-lon_dot*cos(lat))*diff(conj(vd),vd),                                                                                                     -lon_dot*cos(lat)*conj(vn)+lon_dot*sin(lat)*conj(vd),                                                                                                                                                        0,                                                                                                                                                        0]
[                                                                                            cos(phi)*cos(theta)*conj(A_yb)-sin(phi)*cos(theta)*conj(A_zb),                                                                     -cos(theta)*conj(A_xb)-sin(phi)*sin(theta)*conj(A_yb)-cos(phi)*sin(theta)*conj(A_zb),                                                                                                                                                        0,                                                                                                                                lat_dot*diff(conj(vn),vn),                                                                               (5380631166442341/36893488147419103232+lon_dot*cos(lat))*diff(conj(ve),ve),                                                                                                                                                        0,                                                                                                                               -lon_dot*sin(lat)*conj(ve),                                                                                                                                                        0,                                                                                                                                                        0]
[                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                        1/(Rmeridian+hgt),                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                    -vn/(Rmeridian+hgt)^2]
[                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                 1/(Rnormal+hgt)/cos(lat),                                                                                                                                                        0,                                                                                                                     ve/(Rnormal+hgt)/cos(lat)^2*sin(lat),                                                                                                                                                        0,                                                                                                                             -ve/(Rnormal+hgt)^2/cos(lat)]
[                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        1,                                                                                                                                                        0,                                                                                                                                                        0,                                                                                                                                                        0]




This one is in terms of V ned (ie assuming these already calculated by the INS)