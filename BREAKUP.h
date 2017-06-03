//************************** BREAKUP ***************************
// Determines the breakup kinematics (energies(2,3), momenta(0,1), angles(4)), 
// given the (i) breakup angle of light ion, 
//          (ii) total energy of two ions, and
//         (iii) momentum of initial tgt nucleus A.
// NOTE: angl1 is ANTI-clockwize from +x-axis!
//       angl2 is CLOCKwize from +x-axis!
//

float *breakup(Float_t p_0, Float_t m_0, Float_t m_1, Float_t m_2, Float_t E_0, Float_t E_st, Float_t Q);
float MAX_CONE(float m1, float m2, float Et, float p0);

