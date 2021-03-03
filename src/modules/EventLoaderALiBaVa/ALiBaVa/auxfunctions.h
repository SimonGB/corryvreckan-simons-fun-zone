#include "DataFileRoot.h"
#include <TH1.h>

TH1 *draw_gain_hist(DataFileRoot &A,
                      const char *hname, double factor, bool flip, int chip, int odd, int mxevts);

void save_text_file(TH1 *h1, const char *name);

int ALiBaVa_loader(DataFileRoot *A,
                       const char *data_file, const char *cal_file, const char *ped_file);
