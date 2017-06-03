 // Two-body Kinematics code to calculate: 
 // outgoing energy 
 // based on input: mass of projectile, target, ejectile and residual nucleus,  
 // incident energy, Q-value, excitation energy and scattering angle [in radians].
 // 
 //                                              mb
 //                                            /
 // reaction is: a + A -> B + b    --ma->--- mA
 //                                            \\
 //                                             \\     
 //                                              mB
 // JJvZ - Jan 2013
 //

float E_OUT(Float_t Ein, Float_t mp, Float_t mtgt, Float_t me, Float_t mB, Float_t Qgs, Float_t Estar, Float_t THscat);

