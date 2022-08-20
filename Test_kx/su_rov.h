#ifndef SU_ROV_H
#define SU_ROV_H

#include <QObject>
#include <QVector>
#include "qkx_coeffs.h"
#include "kx_protocol.h"
#include <QTimer>

#define ANPA_MOD_CNT 24

extern double X[2000][2];
extern QVector<double> K;

struct InitData {
    double k_gamma;
    double m;
    double Farx[3];
    double delta_m;
    double cv1[4];
    double cv2[4];
    double cw1[4];
    double cw2[4];
    double lambda[6][6];
    double Ta[6][8];
    double J[3];
    double kd;
    double hx;  //metacentricheskaya vysota
    double hy;  //metacentricheskaya vysota
    double hz;  //metacentricheskaya vysota
    double Td;
    double depth_limit;
    double max_depth;
}; //struct InitData

class SU_ROV : public QObject {
    Q_OBJECT
public:
    explicit SU_ROV(QObject *parent = 0);
    virtual ~SU_ROV();
signals:

public slots:
private:
    void start();

public:
    void model(const float Upnp,const float Upnl,const float Uznp,const float Uznl,
               const float Upvp, const float Upvl, const float Uzvl, const float Uzvp);
    void runge(const float Upnp,const float Upnl,const float Uznp,const float Uznl,
               const float Upvp, const float Upvl, const float Uzvl, const float Uzvp,const float Ttimer,const float dt=0.01);

    double a[ANPA_MOD_CNT];
    double da[ANPA_MOD_CNT];
    double delta_f;
    //константы
    double k_gamma;
    double m;
    double Farx[3];
    double g;
    double G;
    double delta_m;
    double cv1[4];
    double cv2[4];
    double cw1[4];
    double cw2[4];
    double lambda[7][7];
    double Ta[7][9];
    double C[7][7];
    double Vt[7];
    double Wv[7]; //вектор силы моментов, вызванных внешними возмущениями
    double J[4];
    double kd;
    double hx;   //metacentricheskaya vysota
    double hy;   //metacentricheskaya vysota
    double hz;   //metacentricheskaya vysota
    double Td;
    double depth_limit;
    double max_depth;
    //переменные
    double sumX, sumZ;
    double cur_depth, Wx, Wy, Wz;
    double Psi_g, Gamma_g, Tetta_g;

    double Psi_gi, W_Psi_g, W_Gamma_g, W_Tetta_g;
    int N;
    double deltaSx, deltaSz;

    double Plz,Plp,Pmpr,Pml; //упоры движителей
    double Ppnp, Ppnl, Pznp, Pznl, Ppvp, Ppvl, Pzvl, Pzvp;
    double Ppnp_x, Ppnl_x, Pznp_x, Pznl_x, Ppvp_x, Ppvl_x, Pzvl_x, Pzvp_x;
    double Ppnp_y, Ppnl_y, Pznp_y, Pznl_y, Ppvp_y, Ppvl_y, Pzvl_y, Pzvp_y;
    double Ppnp_z, Ppnl_z, Pznp_z, Pznl_z, Ppvp_z, Ppvl_z, Pzvl_z, Pzvp_z;
    double Upnp, Upnl, Uznp, Uznl, Upvp, Upvl, Uzvl, Uzvp; //напряжения движителей

    double FloatageX, FloatageY, FloatageZ, Fdx, Fdy, Fdz, Fgx, Fgy, Fgz, Fcx, Fcy, Fcz;
    double Mdx, Mdy, Mdz, Mgx, Mgy, Mgz, Mcx, Mcy, Mcz;
    double Mpnp_x, Mpnl_x, Mznp_x, Mznl_x, Mpvp_x, Mpvl_x, Mzvl_x, Mzvp_x;
    double Mpnp_y, Mpnl_y, Mznp_y, Mznl_y, Mpvp_y, Mpvl_y, Mzvl_y, Mzvp_y;
    double Mpnp_z, Mpnl_z, Mznp_z, Mznl_z, Mpvp_z, Mpvl_z, Mzvl_z, Mzvp_z;
    double Max,May,Maz; // моменты от силы Архимеда

    double x_global, y_global, z_global;
    double vx_local,  vy_local, vz_local;  //lineinye skorosti SPA v svyazannyh osyah
    double vx_global, vy_global, vz_global;

public:
    void resetModel();
    void tick(const float Upnp,const float Upnl,const float Uznp,const float Uznl,
              const float Upvp, const float Upvl, const float Uzvl, const float Uzvp,const float Ttimer);
    float Fx,Fy,Fz; //total forces for XYZ-axis
    float Mx,My,Mz; //total moments for XYZ-axis
protected:
    Qkx_coeffs * K_protocol=nullptr;
    x_protocol * X_protocol=nullptr;
    QTimer timer;
};

#endif // SU_ROV_H
