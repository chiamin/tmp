#ifndef __NONINTERCHAIN_H_CMC__
#define __NONINTERCHAIN_H_CMC__
#include "itensor/all.h"
#include "ReadWriteFile.h"
using namespace itensor;
using namespace std;

template <typename T>
using OpTerm1T = tuple <T, string, int>;
template <typename T>
using OpTerm2T = tuple <T, string, int, string, int>;
using OpTerm1 = OpTerm1T<Real>;
using OpTerm2 = OpTerm2T<Real>;
using OpTerm1C = OpTerm1T<Cplx>;
using OpTerm2C = OpTerm2T<Cplx>;

class NonInterChain
{
    public:
        NonInterChain () {}
        NonInterChain (const Matrix& H0, string name)
        : _name (name)
        {
            diagHermitian (H0, _Uik, _ens);
            for(int i = 0; i < _ens.size(); i++)
                _ops.emplace_back (_ens(i),"N",i+1);
        }

        void add_mu (Real mu)
        {
            for(int i = 0; i < _ens.size(); i++)
                _ops.emplace_back (-mu,"N",i+1);
        }

        const string&          name ()      const { return _name; }
        const Matrix&          Uik  ()      const { return _Uik; }
              Vector           Ui   (int i) const { mycheck (i > 0 and i <= nrows(_Uik), "out of range");
                                                    return Vector (row (_Uik, i-1)); }
        const Vector&          ens  ()      const { return _ens; }
        const int              L    ()      const { return _ens.size(); }
        const vector<OpTerm1>& ops  ()      const { return _ops; }

        void write (ostream& s) const
        {
            itensor::write(s,_Uik);
            itensor::write(s,_ens);
            iutility::write(s,_ops);
            iutility::write(s,_name);
        }
        void read (istream& s)
        {
            itensor::read(s,_Uik);
            itensor::read(s,_ens);
            iutility::read(s,_ops);
            iutility::read(s,_name);
        }

    private:
        Matrix _Uik;
        Vector _ens;
        vector<OpTerm1> _ops;
        string _name;
};
#endif
