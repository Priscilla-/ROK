//
//  polarizations.c
//  X_ROK_V1.7
//
//  Created by Priscilla on 21/07/2013.
//
//

#include <stdio.h>

struct observer                                       // Structure that contains all the information related to the Observers (which are assumed to be located far away from the Source)
{                                                     //
	double d_l;                                         // Luminosity Distance from the Source to the Observer
    //
	double n_sbs[4];                                    // Normal vector pointing from the Source to the Observer in the Solar Barycentric System
	double n_src[4];                                    // Normal vector pointing from the Source to the Observer in the Source Frame
    //
    double r_sbs[4];                                    // Normal vector pointing from the Barycenter of the Solar System to the Source in the Solar Barycentric System
	double s_sbs[4];                                    // Unit Vector that Points in the MBH Spin written in the SBS reference frame
	double p_sbs[4], q_sbs[4];                          // Normal vectors (in the Solar Barycentric System) that form, together with n, a spatial orthonormal triad with respect to the Euclidean scalar product defined via the Kronecker delta (independently of the Cartesian coordinates choosen)
    //
    double s_src[4];                                    // Unit Vector that Points in the MBH Spin written in the Source reference frame
    double p_src[4], q_src[4];                          // Normal vectors (in the Source Reference Frame) that form, together with n, a spatial orthonormal triad with respect to the Euclidean scalar product defined via the Kronecker delta (independently of the Cartesian coordinates choosen)
    //
    double polarization_tensor_plus[4][4];              // Polarization Tensor +
    double polarization_tensor_cross[4][4];             // Polarization Tensor x
    //
	double h_plus;                                      // Multipolar Waveform + polarization associated with the Observer
	double h_cross;                                     // Multipolar Waveform x polarization associated with the Observer
    //
    double LISAtrig[10];                                // Some constant quantities that facilitate the computation of the LISA Response Function
};

struct observer GV_LISA;                              // Observatory



/*----- Setting Time Independent Quantities for GW observations -----*/
/**----- Distance from the Solar System Barycenter to the Source/EMRI (D_L) in MBH mass units or D_L/mu-----**/
factor_GParsecs_to_MBHmass = (1.0e9*PARSEC)/GC_Length_conversion_factor;

if (GV_EMRI_Offset_of_Parameter[3] > 1.0e-12)
GV_LISA.d_l = GV_D_L*factor_GParsecs_to_MBHmass/GV_Mass_SCO;
else
GV_LISA.d_l = Dl_d_SCO*factor_GParsecs_to_MBHmass;

/**----- Unit Vector that Points from the SBS to the Source written in the SBS reference frame -----**/
GV_LISA.r_sbs[1] = sin(GV_Theta_Source)*cos(GV_Phi_Source);
GV_LISA.r_sbs[2] = sin(GV_Theta_Source)*sin(GV_Phi_Source);
GV_LISA.r_sbs[3] = cos(GV_Theta_Source);

/**----- Unit Vector that Points in the MBH Spin written in the SBS reference frame -----**/
GV_LISA.s_sbs[1] = sin(GV_Theta_Kerr)*cos(GV_Phi_Kerr);
GV_LISA.s_sbs[2] = sin(GV_Theta_Kerr)*sin(GV_Phi_Kerr);
GV_LISA.s_sbs[3] = cos(GV_Theta_Kerr);

/**----- Unit Vector that Points in the MBH Spin written in the Source reference frame -----**/
GV_LISA.s_src[1] = 0.0;
GV_LISA.s_src[2] = 0.0;
GV_LISA.s_src[3] = 1.0;

/**----- Unit Vector that points from LISA to the EMRI written in the SBS reference frame [NOTE: This is an approximation.  In general this is a time dependent vector] -----**/
GV_LISA.n_sbs[1] = GV_LISA.r_sbs[1];
GV_LISA.n_sbs[2] = GV_LISA.r_sbs[2];
GV_LISA.n_sbs[3] = GV_LISA.r_sbs[3];

/**----- Unit Vector that points from LISA to the EMRI written in the Source reference frame [NOTE: This is an approximation.  In general this is a time dependent vector] -----**/
GV_LISA.n_src[1] = cos(GV_Theta_Kerr)*cos(GV_Phi_Kerr)*GV_LISA.n_sbs[1] + cos(GV_Theta_Kerr)*sin(GV_Phi_Kerr)*GV_LISA.n_sbs[2] - sin(GV_Theta_Kerr)*GV_LISA.n_sbs[3];
GV_LISA.n_src[2] = -sin(GV_Phi_Kerr)*GV_LISA.n_sbs[1] + cos(GV_Phi_Kerr)*GV_LISA.n_sbs[2];
GV_LISA.n_src[3] = sin(GV_Theta_Kerr)*cos(GV_Phi_Kerr)*GV_LISA.n_sbs[1] + sin(GV_Theta_Kerr)*sin(GV_Phi_Kerr)*GV_LISA.n_sbs[2] + cos(GV_Theta_Kerr)*GV_LISA.n_sbs[3];


/*----- Computing the unit vectors p and q -----*/
/**----- n x s -----**/
vector_product(GV_LISA.p_src, GV_LISA.n_src, GV_LISA.s_src);

/**----- | n x s | -----**/
vector_norm(&norm_of_vector, GV_LISA.p_src);

/**----- p = ( n x s ) / | n x s | -----**/
for (i=1;i<=3;i++)
GV_LISA.p_src[i] /= norm_of_vector;

/**----- q = p x n -----**/
vector_product(GV_LISA.q_src, GV_LISA.p_src, GV_LISA.n_src);

/**----- Computing the Waveform Polarizations in the Source Reference System -----**/
for (i=1;i<=3;i++)
{
    for (j=1;j<=3;j++)
    {
        GV_LISA.polarization_tensor_plus[i][j] = (GV_LISA.p_src[i])*(GV_LISA.p_src[j]) - (GV_LISA.q_src[i])*(GV_LISA.q_src[j]);
        GV_LISA.polarization_tensor_cross[i][j] = (GV_LISA.p_src[i])*(GV_LISA.q_src[j]) + (GV_LISA.q_src[i])*(GV_LISA.p_src[j]);
    }
}

