#include "rov_model.h"


ROV_Model::ROV_Model(QObject *parent) : QObject(parent) {
    resetModel();
    k_gamma = 0.3;
    m = 100;
    delta_m = 2;
    cv1[1] = 109; cv1[2] = 950; cv1[3] = 633;
    cv2[1] = 10.9; cv2[2] = 114; cv2[3] = 76;
    cw1[1] = 228.6; cw1[2] = 366; cw1[3] = 366; // kak v rabote Egorova
    cw2[1] = 2.29; cw2[2] = 36.6; cw2[3] = 36.6;
    Vt[1] = 1;Vt[2] = 1; Vt[3] = 1; Vt[4] = 0; Vt[5] = 0; Vt[6] = 0;
    lambda[1][1] = 50; lambda[2][2] = 101; lambda[3][3] = 101;
    lambda[4][4] = 50; lambda[5][5] = 50; lambda[6][6] = 50;
    Ta[1][1] = 0.7071; Ta[1][2] = 0.7071; Ta[1][3] = -0.7071; Ta[1][4] = -0.7071; Ta[1][5] = 0.7071; Ta[1][6] = 0.7071; Ta[1][7] = -0.7071; Ta[1][8] = -0.7071;
    Ta[2][1] = -0.7071; Ta[2][2] = 0.7071; Ta[2][3] = 0.7071; Ta[2][4] = -0.7071; Ta[2][5] = -0.7071; Ta[2][6] = 0.7071; Ta[2][7] = 0.7071; Ta[2][8] = -0.7071;
    Ta[3][1] = 0.7071; Ta[3][2] = 0.7071; Ta[3][3] = 0.7071; Ta[3][4] = 0.7071; Ta[3][5] = -0.7071; Ta[3][6] = -0.7071; Ta[3][7] = -0.7071; Ta[3][8] = -0.7071;
    Ta[4][1] = 117.5; Ta[4][2] = -117.5; Ta[4][3] = -117.5; Ta[4][4] = 117.5; Ta[4][5] = -117.5; Ta[4][6] = 117.5; Ta[4][7] = 117.5; Ta[4][8] = -117.5;
    Ta[5][1] = -3.83; Ta[5][2] = -3.83; Ta[5][3] = 3.83; Ta[5][4] = 3.83; Ta[5][5] = 3.83; Ta[5][6] = 3.83; Ta[5][7] = -3.83; Ta[5][8] = -3.83;
    Ta[6][1] = -121.33; Ta[6][2] = 121.33; Ta[6][3] = -121.33; Ta[6][4] = 121.33; Ta[6][5] = -121.33; Ta[6][6] = 121.33; Ta[6][7] = -121.33; Ta[6][8] = 121.33;
    //матрица сил и моментов инерции (проверить вторую матрицу, пока я вбила из мат модели, но кажется там я ошиблась она же не симметрична относительно оси, что странно)
    C[1][1] = 0; C[1][2] = (m+lambda[2][2])*da[6]; C[1][3] = -(m + lambda[3][3])*da[5]; C[1][4] = 0; C[1][5] = 0; C[1][6] = 0;
    C[2][1] = -(m + lambda[1][1])*da[6]; C[2][2] = 0; C[2][3] = (m + lambda[3][3])*da[4]; C[2][4] = 0; C[2][5] = 0; C[2][6] = 0;
    C[3][1] = (m + lambda[1][1])*da[5]; C[3][2] = -(m+lambda[2][2])*da[4]; C[3][3] = 0; C[3][4] = 0; C[3][5] = 0; C[3][6] = 0;
    C[4][1] = 0; C[4][2] = 0; C[4][3] = 0; C[4][4] = 0; C[4][5] = -(J[3]+lambda[6][6])*da[6]; C[4][6] = (J[2]+lambda[5][5])*da[5];
    C[5][1] = 0; C[5][2] = 0; C[5][3] = 0; C[5][4] = (J[3]+lambda[6][6])*da[6]; C[5][5] = 0; C[5][6] = -(J[1]+lambda[4][4])*da[4];
    C[6][1] = 0; C[6][2] = 0; C[6][3] = 0; C[6][4] = -(J[2]+lambda[5][5])*da[5]; C[6][5] = (J[1]+lambda[4][4])*da[4]; C[6][6] = 0;
    J[1] = 4; J[2] = 19.8; J[3] = 19.8; //moment inercii apparata vdol sootvetstvuushih osei
    kd = 3; //koefficient usilenija dvizhitelei
    h = 0.018; //metacentricheskaya vysota
    Td = 0.15; //postojannaya vremeni dvizhitelei
    //koordinaty uporov dvizhitelei otnositelno centra mass apparata
    l1=0.14;
    l2=0.14;
    depth_limit=100;
    max_depth=100;
}

void ROV_Model::model(const float Umvl,const float Umnl,const float Umvp,const float Umnp) {
    int limit1, limit2;
    double G,delta_f;

    //модули упоров движителей
    Ppnp = a[7];  // передний нижний правый(1)
    Ppnl = a[8];  // передний нижний левый(2)
    Pznp = a[9];  // задний нижний левый(3)
    Pznl = a[10];  //задний нижний правый(4)
    Ppvp = a[11];  // передний верхний правый(5)
    Ppvl = a[12];  // передний верхний левый(6)
    Pzvl = a[13];  // задний верхний левый(7)
    Pzvp = a[14];  //задний верхний правый(8)

    //проекции упоров движителей на продольную ось апарата X
    Ppnp_x = Ppnp*Ta[1][1];
    Ppnl_x = Ppnl*Ta[1][2];
    Pznp_x = Pznp*Ta[1][3];
    Pznl_x = Pznl*Ta[1][4];
    Ppvp_x = Ppvp*Ta[1][5];
    Ppvl_x = Ppvl*Ta[1][6];
    Pzvl_x = Pzvl*Ta[1][7];
    Pzvp_x = Pzvp*Ta[1][8];

    //проекции упоров движителей на продольную ось апарата Y
    Ppnp_y = Ppnp*Ta[2][1];
    Ppnl_y = Ppnl*Ta[2][2];
    Pznp_y = Pznp*Ta[2][3];
    Pznl_y = Pznl*Ta[2][4];
    Ppvp_y = Ppvp*Ta[2][5];
    Ppvl_y = Ppvl*Ta[2][6];
    Pzvl_y = Pzvl*Ta[2][7];
    Pzvp_y = Pzvp*Ta[2][8];

    //проекции упоров движителей на продольную ось апарата Z
    Ppnp_z = Ppnp*Ta[3][1];
    Ppnl_z = Ppnl*Ta[3][2];
    Pznp_z = Pznp*Ta[3][3];
    Pznl_z = Pznl*Ta[3][4];
    Ppvp_z = Ppvp*Ta[3][5];
    Ppvl_z = Ppvl*Ta[3][6];
    Pzvl_z = Pzvl*Ta[3][7];
    Pzvp_z = Pzvp*Ta[3][8];

    //момент создаваемый движетельным комплексом вокруг оси X
    Mpnp_x = Ppnp*Ta[4][1];
    Mpnl_x = Ppnl*Ta[4][2];
    Mznp_x = Pznp*Ta[4][3];
    Mznl_x = Pznl*Ta[4][4];
    Mpvp_x = Ppvp*Ta[4][5];
    Mpvl_x = Ppvl*Ta[4][6];
    Mzvl_x = Pzvl*Ta[4][7];
    Mzvp_x = Pzvp*Ta[4][8];

    //момент создаваемый движетельным комплексом вокруг оси Y
    Mpnp_y = Ppnp*Ta[5][1];
    Mpnl_y = Ppnl*Ta[5][2];
    Mznp_y = Pznp*Ta[5][3];
    Mznl_y = Pznl*Ta[5][4];
    Mpvp_y = Ppvp*Ta[5][5];
    Mpvl_y = Ppvl*Ta[5][6];
    Mzvl_y = Pzvl*Ta[5][7];
    Mzvp_y = Pzvp*Ta[5][8];

    //момент создаваемый движетельным комплексом вокруг оси Z
    Mpnp_z = Ppnp*Ta[6][1];
    Mpnl_z = Ppnl*Ta[6][2];
    Mznp_z = Pznp*Ta[6][3];
    Mznl_z = Pznl*Ta[6][4];
    Mpvp_z = Ppvp*Ta[6][5];
    Mpvl_z = Ppvl*Ta[6][6];
    Mzvl_z = Pzvl*Ta[6][7];
    Mzvp_z = Pzvp*Ta[6][8];

    double g = 9.81;
    G = m*g; //вес аппарата
    delta_f = delta_m * g; //плавучесть (H)
    Farx[1][1] = 0; Farx [1][2] = 0; Farx[1][3] = -10; //коэффициент 10 написан по приколу

    //obnulenie verticalnoi polozhitelnoi skorosti apparata pri dostizhenii poverhnosti
    limit1 = limit2 = 0;
    if (a[15] >= max_depth) {
      a[15] = max_depth;
        if (a[3] <= 0) {
          a[3] = 0;
          limit1 = 1;
      }
    };

    //obnulenie verticalnoi polozhitelnoi skorosti apparata pri dostizhenii dna
    if (a[15] <= 0)
    {
      a[15] = 0;
        if (a[3] >= 0)
      {
          a[3] = 0;
          limit2 = 1;
      }
    };

    Fdx = Ppnp_x + Ppnl_x + Pznp_x + Pznl_x + Ppvp_x + Ppvl_x + Pzvl_x + Pzvp_x; // вектор сил и моментов, создаваемых движительным комплексом
    Fgx = -cv1[1] * a[1] * fabs(a[1]) - cv2[1] * a[1]; //произведение D1*Vx
    FloatageX = -sin(a[6]) * (G - Farx[1][3]);
    Fcx = C[1][1]*a[1] + C[1][2]*a[2]+C[1][3]*a[3]+C[1][4]*da[4]+C[1][5]*da[5] + C[1][6]*da[6];
    //FloatageX = 0; //обнуление остаточной плавучести
    da[1] = (1/(m + lambda[1][1])) * (Fdx - Fgx -Fcx - FloatageX+Vt[1]); //vx'

    Fdy = Ppnp_y + Ppnl_y + Pznp_y + Pznl_y + Ppvp_y + Ppvl_y + Pzvl_y + Pzvp_y; // вектор сил и моментов, создаваемых движительным комплексом
    Fgy = -cv1[2] * a[2] * fabs(a[2]) - cv2[2] * a[2]; //произведение D1*Vy
    FloatageY = cos(a[6]) * sin(a[5]) * (G - Farx[1][3]);
    Fcy = C[2][1]*a[1] + C[2][2]*a[2]+C[2][3]*a[3]+C[2][4]*da[4]+C[2][5]*da[5] + C[2][6]*da[6];
    //FloatageY = 0; //обнуление остаточной плавучести
    da[2] = (1/(m + lambda[2][2])) * (Fdy - Fgy -Fcy - FloatageY+Vt[2]); //vy'

    Fdz = Ppnp_z + Ppnl_z + Pznp_z + Pznl_z + Ppvp_z + Ppvl_z + Pzvl_z + Pzvp_z; // вектор сил и моментов, создаваемых движительным комплексом
    Fgz = -cv1[3] * a[3] * fabs(a[3]) - cv2[3] * a[3]; //произведение D1*Vz
    FloatageZ = cos(a[6]) * cos(a[5]) * (G - Farx[1][3]);
    Fcz = C[3][1]*a[1] + C[3][2]*a[2]+C[3][3]*a[3]+C[3][4]*da[4]+C[3][5]*da[5] + C[3][6]*da[6];
    //FloatageZ = 0; //обнуление остаточной плавучести
    da[3] = (1/(m + lambda[3][3])) * (Fdz - Fgz -Fcz - FloatageZ+Vt[3]); //vz'
// da[4] -> производная угла крена da[5]-> производная угла дифферента da[6]-> производная дкурса,
//следующие 3 уравнения это Кинематические уравнения для углов Эйлера-Крылова
//описывающее преобразование вектора угловых скоростей относительно осей НПА Ox,Oy,Oz , в вектор
//угловых скоростей  по курсу, дифференту и крену соответственно.

    da[4] = a[17] + (1/cos(a[5]) * ((a[18]) * sin(a[4]) * sin(a[5])  - sin(a[5]) * cos(a[4]) * a[19]));  //proizvodnaya krena

    da[5] = ((a[18]) * cos(a[4]) - sin(a[4]) * a[19]);  //proizvodnaya differenta

    da[6] = (1/cos(a[5])) * (a[18] * sin(a[4]) + cos(a[4]) * (a[18])); //proizvodnaya kursa

 //пока не шарю что это {
 //   X[17][0]=da[7] = (1/Td) * (kd * (double)Umvp - Pmvp);  // маршевый верхний правый

  //  da[8] = (1/Td) * (kd * (double)Umvl - Pmvl); //маршевый верхний левый

  //  da[9] = (1/Td) * (kd * (double)Umnp - Pmnp);  // маршевый нижний правый

  //  da[10] = (1/Td) * (kd * (double)Umnl - Pmnl);  //маршевый нижний левый

 //   da[11] = 0;

   // da[12] = 0;

  //  da[13] = 0; }

    double alfa[4][4]; //матрица перевода из связанной СК в глобальную СК
    a[4] = -a[4];
    alfa[1][1] = cos(a[5])*cos(a[6]);
    alfa[2][1] = sin(a[6])*cos(a[5]);
    alfa[3][1] = -sin(a[5]);
    alfa[1][2] = cos(a[6])*sin(a[5])*sin(a[4])-cos(a[4])*sin(a[6]);
    alfa[2][2] = cos(a[6])*cos(a[4])+sin(a[4])*sin(a[5])*sin(a[6]);
    alfa[3][2] = sin(a[4])*cos(a[5]);
    alfa[1][3] = sin(a[6])*sin(a[4])+cos(a[6])*cos(a[4])*sin(a[5]);
    alfa[2][3] = sin(a[5])*sin(a[6])*cos(a[4])-cos(a[6])*sin(a[4]);
    alfa[3][3] = cos(a[5])*cos(a[4]);
    a[4] = -a[4];
//кажется тут нужно думать со смещением, вся модель была сделана максимум на 7 движетелей, поэтому a[14] накладывается друг на друга, если что завтра этим займусь
    da[14] = alfa[1][1] * a[1] + alfa[1][2] * a[2] + alfa[1][3] * a[3];
    //dx_global

    da[15] = alfa[2][1] * a[1] + alfa[2][2] * a[2] + alfa[2][3] * a[3];
    //dy_global

    da[16] = alfa[3][1] * a[1] + alfa[3][2] * a[2] + alfa[3][3] * a[3];
    //dz_global

    double Fa = G + delta_f;
    double Fax = sin(a[6])*Fa;
    //float Fay = cos(a[5])*cos(a[6])*Fa;
    double Faz = -sin(a[5])*cos(a[6])*Fa;

    Mdx = k_gamma*( -Pmvp_x - Pmnl_x + Pmnp_x + Pmvl_x );
    Mgx = -cw1[0] * a[17] * fabs(a[17]) - cw2[0] * a[17];
    Max = Faz*h;
    //Max = 0; //obnulenie momenta ot sily Arhimeda
    da[17] = (1/(J[0] + lambda[4][4])) * (Mdx + Mgx + Max);

    Mdy = l2*(-Pmvp_x + Pmvl_x + Pmnl_x - Pmnp_x);
    Mgy = -cw1[1] * a[18] * fabs(a[18]) - cw2[1] * a[18];
    da[18] = (1/(J[1] + lambda[5][5])) * (Mdy + Mgy);

    Mdz = l1*(-Pmvp_x - Pmvl_x + Pmnl_x + Pmnp_x);
    Mgz = -cw1[2] * a[19] * fabs(a[19]) - cw2[2] * a[19];
    Maz = -h*Fax;
    //Maz = 0; //obnulenie momenta ot sily Arhimeda
    da[19] = (1/(J[2] + lambda[6][6])) * (Mdz + Mgz +Maz);

    da[20] = a[1];
    da[21] = a[2];
    da[22] = a[3];

  

}


void ROV_Model::resetModel(){
    for (int i=0;i<ANPA_MOD_CNT;i++) {a[i] = 0.0f; da[i]=0.0f;}
}

void ROV_Model::tick(const float Umvl,const float Umnl,
                     const float Umvp,const float Umnp,const float Ttimer){
    runge(Umvl,Umnl,Umvp,Umnp,Ttimer,Ttimer);
}

ROV_Model::~ROV_Model(){

}

void ROV_Model::runge(const float Umvl, const float Umnl, const float Umvp, const float Umnp, const float Ttimer, const float dt) {
    const double Kc = 180/M_PI;
    double a1[23], y[23];
    int i;
    const double H1 = dt;
    const int n = ANPA_MOD_CNT;
    model(Umvl,Umnl,Umvp,Umnp);
    for (i = 1; i < n; i++) {
      a1[i] = a[i];
      y[i] = da[i];
      a[i] = a1[i] + 0.5 * H1 * da[i];
    }

    model(Umvl,Umnl,Umvp,Umnp);
    for (i = 1; i < n; i++)
    {
      y[i] = y[i]+ 2 * da[i];
      a[i] = a1[i] + 0.5 * H1 * da[i];
    }
    model(Umvl,Umnl,Umvp,Umnp);
    for (i = 1; i < n; i++) {
      y[i] = y[i] + 2 * da[i];
      a[i] = a1[i] + H1 * da[i];
    }
    model(Umvl,Umnl,Umvp,Umnp);
    for (i = 1; i < n; i++) {
      a[i] = a1[i] + (H1 / 6) * (y[i] + da[i]);
    }


    //данные в СУ ( с преобразованием координат)

    x_global = a[15]; //koordinata apparata v globalnoi SK
    y_global = a[16];  //otstojanie ot dna otnositelno repernoi tochki, kotoraja na dne
    cur_depth = max_depth - y_global;  //tekush"aya glubina SPA
    z_global = a[17]; //koordinaty apparata v globalnoi SK (преобразование координат)
    Wx = a[18] * Kc; //uglovye skorosti SPA v svyazannyh osyah v gradus/sekunda
    Wy = a[19] * Kc;
    Wz = a[20] * Kc;

    vx_local = a[1]; vy_local = a[2]; vz_local = a[3];  //lineinye skorosti SPA v svyazannyh osyah
    vx_global = da[14]; vy_global = da[15]; vz_global = da[16];  // lineinye skorosti SPA v globalnyh osyah

    Gamma_g = a[4] * Kc; // ugol krena
    Tetta_g = a[5] * Kc; // ugol differenta
    Psi_g = a[6] * Kc; // ugol kursa (преобразование координат)

    W_Psi_g = da[4] * Kc; // proizvodnaya ugla kursa
    W_Gamma_g = da[5] * Kc; // proizvodnaya ugla krena
    W_Tetta_g = da[6] * Kc; // proizvodnaya ugla differenta

    N = fabs(Psi_g / 360);
    if (Psi_g >= 360) Psi_gi = Psi_g - N * 360; // ugol kursa na indikaciu
    if (Psi_g <= -360) Psi_gi = Psi_g + N * 360;

    deltaSx = vx_local * Ttimer; //prirash"enie koordinaty X dlya SVS (v svyazannoi s SPA SK)
    sumX += deltaSx;

    deltaSz = vz_local * Ttimer; //prirash"enie koordinaty Z dlya SVS (v svyazannoi s SPA SK)
    sumZ += deltaSz;


}



