#include <QCoreApplication>
#include "su_rov.h"

double X[2000][2];

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    SU_ROV su;
    su.tick(7,10,8,7,0,0,0,0,0.1);
    return a.exec();
}
