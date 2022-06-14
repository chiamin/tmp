#ifndef __SINGLEPARTICLE_H_CMC__
#define __SINGLEPARTICLE_H_CMC__

inline int mix_ind (int xi, int yi, int Lx, int Ly, int first_index=1, int first_x_index=0)
{
    return (xi-first_x_index) * Ly + (yi-first_x_index) + first_index;
}

Matrix Hamilt_k (int L, Real t, Real mu, Real damp_fac=1., bool damp_from_right=true, bool verbose=false)
{
    cout << "L = " << L << endl;
    Matrix H (L,L);
    for(int i = 0; i < L; i++)
    {
        H(i,i) = -mu;
        if (i != L-1)
        {
            int damp_dist = (damp_from_right ? L-2-i : i);
            Real ti = t * pow (damp_fac, damp_dist);
            H(i,i+1) = -ti;
            H(i+1,i) = -ti;
            if (verbose)
                cout << "Hk, t " << i << " = " << ti << endl;
        }
    }
    return H;
}

#endif
