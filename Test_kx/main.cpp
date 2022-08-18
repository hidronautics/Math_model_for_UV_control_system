#include <QCoreApplication>
#include "su_rov.h"
#include "rov_model.h"

double X[2000][2];

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    SU_ROV su;
    ROV_Model rov;
    return a.exec();
}
