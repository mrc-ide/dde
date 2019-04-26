#include "dopri_853.h"

// *Massive* constant list
#define  C2      0.526001519587677318785587544488e-01
#define  C3      0.789002279381515978178381316732e-01
#define  C4      0.118350341907227396726757197510
#define  C5      0.281649658092772603273242802490
#define  C6      0.333333333333333333333333333333
#define  C7      0.25
#define  C8      0.307692307692307692307692307692
#define  C9      0.651282051282051282051282051282
#define  C10     0.6
#define  C11     0.857142857142857142857142857142
#define  C14     0.1
#define  C15     0.2
#define  C16     0.777777777777777777777777777778
#define  B1      5.42937341165687622380535766363e-02
#define  B6      4.45031289275240888144113950566
#define  B7      1.89151789931450038304281599044
#define  B8     -5.8012039600105847814672114227
#define  B9      3.1116436695781989440891606237e-01
#define  B10    -1.52160949662516078556178806805e-01
#define  B11     2.01365400804030348374776537501e-01
#define  B12     4.47106157277725905176885569043e-02
#define  BHH1    0.244094488188976377952755905512
#define  BHH2    0.733846688281611857341361741547
#define  BHH3    0.220588235294117647058823529412e-01
#define  ER1     0.1312004499419488073250102996e-01
#define  ER6    -0.1225156446376204440720569753e+01
#define  ER7    -0.4957589496572501915214079952
#define  ER8     0.1664377182454986536961530415e+01
#define  ER9    -0.3503288487499736816886487290
#define  ER10    0.3341791187130174790297318841
#define  ER11    0.8192320648511571246570742613e-01
#define  ER12   -0.2235530786388629525884427845e-01
#define  A21     5.26001519587677318785587544488e-02
#define  A31     1.97250569845378994544595329183e-02
#define  A32     5.91751709536136983633785987549e-02
#define  A41     2.95875854768068491816892993775e-02
#define  A43     8.87627564304205475450678981324e-02
#define  A51     2.41365134159266685502369798665e-01
#define  A53    -8.84549479328286085344864962717e-01
#define  A54     9.24834003261792003115737966543e-01
#define  A61     3.7037037037037037037037037037e-02
#define  A64     1.70828608729473871279604482173e-01
#define  A65     1.25467687566822425016691814123e-01
#define  A71     3.7109375e-02
#define  A74     1.70252211019544039314978060272e-01
#define  A75     6.02165389804559606850219397283e-02
#define  A76    -1.7578125e-02
#define  A81     3.70920001185047927108779319836e-02
#define  A84     1.70383925712239993810214054705e-01
#define  A85     1.07262030446373284651809199168e-01
#define  A86    -1.53194377486244017527936158236e-02
#define  A87     8.27378916381402288758473766002e-03
#define  A91     6.24110958716075717114429577812e-01
#define  A94    -3.36089262944694129406857109825
#define  A95    -8.68219346841726006818189891453e-01
#define  A96     2.75920996994467083049415600797e+01
#define  A97     2.01540675504778934086186788979e+01
#define  A98    -4.34898841810699588477366255144e+01
#define  A101    4.77662536438264365890433908527e-01
#define  A104   -2.48811461997166764192642586468
#define  A105   -5.90290826836842996371446475743e-01
#define  A106    2.12300514481811942347288949897e+01
#define  A107    1.52792336328824235832596922938e+01
#define  A108   -3.32882109689848629194453265587e+01
#define  A109   -2.03312017085086261358222928593e-02
#define  A111   -9.3714243008598732571704021658e-01
#define  A114    5.18637242884406370830023853209
#define  A115    1.09143734899672957818500254654
#define  A116   -8.14978701074692612513997267357
#define  A117   -1.85200656599969598641566180701e+01
#define  A118    2.27394870993505042818970056734e+01
#define  A119    2.49360555267965238987089396762
#define  A1110  -3.0467644718982195003823669022
#define  A121    2.27331014751653820792359768449
#define  A124   -1.05344954667372501984066689879e+01
#define  A125   -2.00087205822486249909675718444
#define  A126   -1.79589318631187989172765950534e+01
#define  A127    2.79488845294199600508499808837e+01
#define  A128   -2.85899827713502369474065508674
#define  A129   -8.87285693353062954433549289258
#define  A1210   1.23605671757943030647266201528e+01
#define  A1211   6.43392746015763530355970484046e-01
#define  A141    5.61675022830479523392909219681e-02
#define  A147    2.53500210216624811088794765333e-01
#define  A148   -2.46239037470802489917441475441e-01
#define  A149   -1.24191423263816360469010140626e-01
#define  A1410   1.5329179827876569731206322685e-01
#define  A1411   8.20105229563468988491666602057e-03
#define  A1412   7.56789766054569976138603589584e-03
#define  A1413  -8.298e-03
#define  A151    3.18346481635021405060768473261e-02
#define  A156    2.83009096723667755288322961402e-02
#define  A157    5.35419883074385676223797384372e-02
#define  A158   -5.49237485713909884646569340306e-02
#define  A1511  -1.08347328697249322858509316994e-04
#define  A1512   3.82571090835658412954920192323e-04
#define  A1513  -3.40465008687404560802977114492e-04
#define  A1514   1.41312443674632500278074618366e-01
#define  A161   -4.28896301583791923408573538692e-01
#define  A166   -4.69762141536116384314449447206
#define  A167    7.68342119606259904184240953878
#define  A168    4.06898981839711007970213554331
#define  A169    3.56727187455281109270669543021e-01
#define  A1613  -1.39902416515901462129418009734e-03
#define  A1614   2.9475147891527723389556272149
#define  A1615  -9.15095847217987001081870187138
#define  D41    -0.84289382761090128651353491142e+01
#define  D46     0.56671495351937776962531783590
#define  D47    -0.30689499459498916912797304727e+01
#define  D48     0.23846676565120698287728149680e+01
#define  D49     0.21170345824450282767155149946e+01
#define  D410   -0.87139158377797299206789907490
#define  D411    0.22404374302607882758541771650e+01
#define  D412    0.63157877876946881815570249290
#define  D413   -0.88990336451333310820698117400e-01
#define  D414    0.18148505520854727256656404962e+02
#define  D415   -0.91946323924783554000451984436e+01
#define  D416   -0.44360363875948939664310572000e+01
#define  D51     0.10427508642579134603413151009e+02
#define  D56     0.24228349177525818288430175319e+03
#define  D57     0.16520045171727028198505394887e+03
#define  D58    -0.37454675472269020279518312152e+03
#define  D59    -0.22113666853125306036270938578e+02
#define  D510    0.77334326684722638389603898808e+01
#define  D511   -0.30674084731089398182061213626e+02
#define  D512   -0.93321305264302278729567221706e+01
#define  D513    0.15697238121770843886131091075e+02
#define  D514   -0.31139403219565177677282850411e+02
#define  D515   -0.93529243588444783865713862664e+01
#define  D516    0.35816841486394083752465898540e+02
#define  D61     0.19985053242002433820987653617e+02
#define  D66    -0.38703730874935176555105901742e+03
#define  D67    -0.18917813819516756882830838328e+03
#define  D68     0.52780815920542364900561016686e+03
#define  D69    -0.11573902539959630126141871134e+02
#define  D610    0.68812326946963000169666922661e+01
#define  D611   -0.10006050966910838403183860980e+01
#define  D612    0.77771377980534432092869265740
#define  D613   -0.27782057523535084065932004339e+01
#define  D614   -0.60196695231264120758267380846e+02
#define  D615    0.84320405506677161018159903784e+02
#define  D616    0.11992291136182789328035130030e+02
#define  D71    -0.25693933462703749003312586129e+02
#define  D76    -0.15418974869023643374053993627e+03
#define  D77    -0.23152937917604549567536039109e+03
#define  D78     0.35763911791061412378285349910e+03
#define  D79     0.93405324183624310003907691704e+02
#define  D710   -0.37458323136451633156875139351e+02
#define  D711    0.10409964950896230045147246184e+03
#define  D712    0.29840293426660503123344363579e+02
#define  D713   -0.43533456590011143754432175058e+02
#define  D714    0.96324553959188282948394950600e+02
#define  D715   -0.39177261675615439165231486172e+02
#define  D716   -0.14972683625798562581422125276e+03

void dopri853_step(dopri_data *obj, double h) {
  const double t = obj->t;
  const size_t n = obj->n;
  double
    *k1 = obj->k[0],
    *k2 = obj->k[1],
    *k3 = obj->k[2],
    *k4 = obj->k[3],
    *k5 = obj->k[4],
    *k6 = obj->k[5],
    *k7 = obj->k[6],
    *k8 = obj->k[7],
    *k9 = obj->k[8],
    *k10 = obj->k[9];
  double *y = obj->y, *y1 = obj->y1;

  // TODO: do we need to call target here with y, k1?  Looks like
  // that's probably taken care of for us, as a similar call exists in
  // dopri5.
  for (size_t i = 0; i < n; ++i) { // 22
    y1[i] = y[i] + h * A21 * k1[i];
  }
  dopri_eval(obj, t + C2 * h, y1, k2);

  for (size_t i = 0; i < n; ++i) { // 23
    y1[i] = y[i] + h * (A31 * k1[i] + A32 * k2[i]);
  }
  dopri_eval(obj, t + C3 * h, y1, k3);

  for (size_t i = 0; i < n; ++i) { // 24
    y1[i] = y[i] + h * (A41 * k1[i] + A43 * k3[i]);
  }
  dopri_eval(obj, t + C4 * h, y1, k4);

  for (size_t i = 0; i < n; ++i) { // 25
    y1[i] = y[i] + h * (A51 * k1[i] + A53 * k3[i] + A54 * k4[i]);
  }
  dopri_eval(obj, t + C5 * h, y1, k5);

  for (size_t i = 0; i < n; ++i) { // 26
    y1[i] = y[i] + h * (A61 * k1[i] + A64 * k4[i] + A65 * k5[i]);
  }
  dopri_eval(obj, t + C6 * h, y1, k6);

  for (size_t i = 0; i < n; ++i) { // 27
    y1[i] = y[i] + h * (A71 * k1[i] + A74 * k4[i] + A75 * k5[i] + A76 * k6[i]);
  }
  dopri_eval(obj, t + C7 * h, y1, k7);

  for (size_t i = 0; i < n; ++i) { // 28
    y1[i] = y[i] + h * (A81 * k1[i] + A84 * k4[i] + A85 * k5[i] +
                        A86 * k6[i] + A87 * k7[i]);
  }
  dopri_eval(obj, t + C8 * h, y1, k8);

  for (size_t i = 0; i < n; ++i) { // 29
    y1[i] = y[i] + h * (A91 * k1[i] + A94 * k4[i] + A95 * k5[i] +
                        A96 * k6[i] + A97 * k7[i] + A98 * k8[i]);
  }
  dopri_eval(obj, t + C9 * h, y1, k9);

  for (size_t i = 0; i < n; ++i) { // 30
    y1[i] = y[i] + h * (A101 * k1[i] + A104 * k4[i] + A105 * k5[i] +
                        A106 * k6[i] + A107 * k7[i] + A108 * k8[i] +
                        A109 * k9[i]);
  }
  dopri_eval(obj, t + C10 * h, y1, k10);

  for (size_t i = 0; i < n; ++i) { // 31
    y1[i] = y[i] + h * (A111 * k1[i] + A114 * k4[i] + A115 * k5[i] +
                        A116 * k6[i] + A117 * k7[i] + A118 * k8[i] +
                        A119 * k9[i] + A1110 * k10[i]);
  }
  dopri_eval(obj, t + C11 * h, y1, k2);

  double t_next = t + h;
  for (size_t i = 0; i < n; ++i) { // 32
    y1[i] = y[i] + h * (A121 * k1[i] + A124  * k4[i]  + A125  * k5[i] +
                        A126 * k6[i] + A127  * k7[i]  + A128  * k8[i] +
                        A129 * k9[i] + A1210 * k10[i] + A1211 * k2[i]);
  }
  dopri_eval(obj, t_next, y1, k3);

  for (size_t i = 0; i < n; ++i) { // 35
    k4[i] = B1 * k1[i] + B6 * k6[i] + B7 * k7[i] + B8 * k8[i] +
      B9 * k9[i] + B10 * k10[i] + B11 * k2[i] + B12 * k3[i];
    k5[i] = y[i] + h * k4[i];
  }
}

double dopri853_error(dopri_data *obj) {
  double
    *k1 = obj->k[0],
    *k2 = obj->k[1],
    *k3 = obj->k[2],
    *k4 = obj->k[3],
    *k5 = obj->k[4],
    *k6 = obj->k[5],
    *k7 = obj->k[6],
    *k8 = obj->k[7],
    *k9 = obj->k[8],
    *k10 = obj->k[9];
  double err = 0.0, err2 = 0.0;
  for (size_t i = 0; i < obj->n; ++i) {
    double sk = obj->atol + obj->rtol * fmax(fabs(obj->y[i]), fabs(k5[i]));
    double erri = (k4[i] -
                   BHH1 * k1[i] -
                   BHH2 * k9[i] -
                   BHH3 * k3[i]);
    err2 += square(erri / sk);
    erri = (ER1 * k1[i] +
            ER6 * k6[i] +
            ER7 * k7[i] +
            ER8 * k8[i] +
            ER9 * k9[i] +
            ER10 * k10[i] +
            ER11 * k2[i] +
            ER12 * k3[i]);
    err += square(erri / sk);
  }
  double deno = err + 0.01 * err2;
  // This is some sort of safety catch; it is never triggered in the tests
  deno = deno < 0 ? 1.0 : deno;
  return obj->sign * err * sqrt(1.0 / (obj->n * deno));
}

void dopri853_save_history(dopri_data *obj, double h) {
  double *history = (double*) obj->history->head;
  double *y = obj->y, *y1 = obj->y1;
  double
    *k1 = obj->k[0],
    *k2 = obj->k[1],
    *k3 = obj->k[2],
    *k4 = obj->k[3],
    *k5 = obj->k[4],
    *k6 = obj->k[5],
    *k7 = obj->k[6],
    *k8 = obj->k[7],
    *k9 = obj->k[8],
    *k10 = obj->k[9];
  double t = obj->t;
  size_t n = obj->n;

  // NOTE: We have a function call here, in contrast with dopri5.
  // This call comes from dop853.f:673
  dopri_eval(obj, obj->t + h, k5, k4);

  for (size_t i = 0; i < n; ++i) {
    double ydiff = k5[i] - y[i];
    double bspl = h * k1[i] - ydiff;
    history[             i] = y[i];
    history[    n + i] = ydiff;
    history[2 * n + i] = bspl;
    history[3 * n + i] = ydiff - h * k4[i] - bspl;
    // Next ones are more different than the dopri5 case and
    // significantly uglier:
    history[4 * n + i] = (D41  * k1[i] + D46  * k6[i] + D47  * k7[i]  +
                          D48  * k8[i] + D49  * k9[i] + D410 * k10[i] +
                          D411 * k2[i] + D412 * k3[i]);
    history[5 * n + i] = (D51  * k1[i] + D56  * k6[i] + D57  * k7[i]  +
                          D58  * k8[i] + D59  * k9[i] + D510 * k10[i] +
                          D511 * k2[i] + D512 * k3[i]);
    history[6 * n + i] = (D61  * k1[i] + D66  * k6[i] + D67  * k7[i]  +
                          D68  * k8[i] + D69  * k9[i] + D610 * k10[i] +
                          D611 * k2[i] + D612 * k3[i]);
    history[7 * n + i] = (D71  * k1[i] + D76  * k6[i] + D77  * k7[i]  +
                          D78  * k8[i] + D79  * k9[i] + D710 * k10[i] +
                          D711 * k2[i] + D712 * k3[i]);
  }

  // Then three more function evaluations
  for (size_t i = 0; i < n; ++i) { // 51
    y1[i]=y[i] + h * (A141  * k1[i] + A147  * k7[i]  + A148  * k8[i] +
                      A149  * k9[i] + A1410 * k10[i] + A1411 * k2[i] +
                      A1412 * k3[i] + A1413 * k4[i]);
  }
  dopri_eval(obj, t + C14 * h, y1, k10);
  for (size_t i = 0; i < n; ++i) { // 52
    y1[i]=y[i] + h * (A151  * k1[i] + A156  * k6[i] + A157  * k7[i] +
                      A158  * k8[i] + A1511 * k2[i] + A1512 * k3[i] +
                      A1513 * k4[i] + A1514 * k10[i]);
  }
  dopri_eval(obj, t + C15 * h, y1, k2);
  for (size_t i = 0; i < n; ++i) { // 53
    y1[i]=y[i] + h * (A161  * k1[i]  + A166  * k6[i] + A167  * k7[i] +
                      A168  * k8[i]  + A169  * k9[i] + A1613 * k4[i] +
                      A1614 * k10[i] + A1615 * k2[i]);
  }
  dopri_eval(obj, t + C16 * h, y1, k3);

  // Final history preparation
  for (size_t i = 0; i < n; ++i) {
    history[4 * n + i] = h * (history[4 * n + i] +
                              D413 * k4[i] + D414 * k10[i] +
                              D415 * k2[i] + D416 * k3[i]);
    history[5 * n + i] = h * (history[5 * n + i] +
                              D513 * k4[i] + D514 * k10[i] +
                              D515 * k2[i] + D516 * k3[i]);
    history[6 * n + i] = h * (history[6 * n + i] +
                              D613 * k4[i] + D614 * k10[i] +
                              D615 * k2[i] + D616 * k3[i]);
    history[7 * n + i] = h * (history[7 * n + i] +
                              D713 * k4[i] + D714 * k10[i] +
                              D715 * k2[i] + D716 * k3[i]);
  }

  history[obj->history_idx_time    ] = t;
  history[obj->history_idx_time + 1] = h;
}

double dopri853_interpolate(size_t n, double theta, double theta1,
                            const double *history) {
  double tmp = history[4 * n] + theta *
    (history[5 * n] + theta1 *
     (history[6 * n] + theta *
      history[7 * n]));
  return history[0] + theta *
    (history[n] + theta1 *
     (history[2 * n] + theta *
      (history[3 * n] + theta1 *
       tmp)));
}

bool dopri853_test_stiff(dopri_data *obj, double h) {
  double
    *k3 = obj->k[2],
    *k4 = obj->k[3],
    *k5 = obj->k[4];
  double *y1 = obj->y1;
  double stnum = 0, stden = 0;
  for (size_t i = 0; i < obj->n; ++i) {
    stnum += square(k4[i] - k3[i]);
    stden += square(k5[i] - y1[i]);
  }
  return stden > 0 && fabs(h) * sqrt(stnum / stden) > 6.1;
}
