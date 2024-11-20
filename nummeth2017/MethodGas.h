#pragma once

#include "method.h"


class MethodGas : public Method
{
public:
    virtual void convertToParam(int, Param&);
    virtual void init(); //��������� �������
    virtual void run(); // ������ ����

    ~MethodGas();
protected:

    void initValues(); //���������� � init()

    /*
    pl - ��������� � ����
    pr - ��������� � ����� 

    �� ������ ����������� ���������� fr, fu, fv, fe
                fr = ro * u
    F1 =        fu = ro * u^2 + p
                fv = ro * u * v
                fe = u * (ro * (eps + (u^2 + v^2)/2))


                fr = ro * v
    F2 =        fu = ro * u * v
                fv = ro * v^2 + p
                fe = v * (ro * (eps + (u^2 + v^2)/2))

    ��� ������� ����� �������� �� �������

    */

    void flux(Param pl, Param pr, Vector n, double& fr, double& fu, double& fv, double& fe); //������ �������
    void flux_visc(Param pl, Param pr, Vector n, double& fr, double& fu, double& fv, double& fe, double grad_u_x, double grad_u_y, double grad_v_x, double grad_v_y);


    void bnd(Edge* e, Param p1, Param& p2); // ������������ ����� ��������� ����� (��� ���������)

    double* ro, * ru, * rv, * re;
    double* int_ro, * int_ru, * int_rv, * int_re;

    double* grad_u_x, *grad_u_y, *grad_v_x, *grad_v_y; //-------------------------

    double TMAX; // ����� �� �������� �������
    double TAU; // ��� �� �������

    const double GAM = 1.4;
};


