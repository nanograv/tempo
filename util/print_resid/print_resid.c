/* print_resid.c 
 *
 * Read a resid2.tmp file and do stuff with it
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include "tempo_output.h"

#define MAX_OUTS 10

void usage() {
    printf(
            "Usage: print_resid [output_options] resid_file\n"
            "Options determine what is printed, in what order (max %d):\n"
            "  -m, --mjd        MJD of barycentric TOA\n"
            "  -p, --res_phase  Residual in phase (turns)\n"
            "  -t, --res_sec    Residual in time (sec)\n"
            "  -r, --res_us     Residual in time (us)\n"
            "  -o, --ophase     Orbital phase (turns)\n"
            "  -f, --rf         Barycentric RF (MHz)\n"
            "  -w, --weight     Weight of point in fit\n"
            "  -e, --err        TOA uncertainty (us)\n"
            "  -i, --prefit_sec Pre-fit residual (sec)\n"
            "  -d, --ddm        Delta DM used (pc/cc)\n"
            "Other options:\n"
            "  -s, --stats      Print stats at beginning\n"
            "  -b, --bands      Print blank lines between bands\n"
            "  -h, --help       Print this message\n"
            "Calling with no args is equivalent to:\n"
            "  print_resid -mfreo resid2.tmp\n", 
            MAX_OUTS
          );
}

/* Calc weighted or not rms */
double resid_rms(int npts, struct residual *r, int weight) {
    int i;
    double rsum=0.0, wsum=0.0, w;
    for (i=0; i<npts; i++) {
        w = weight ? r[i].weight : 1.0; 
        rsum += w*r[i].res_sec*r[i].res_sec;
        wsum += w;
    }
    rsum /= wsum;
    rsum = sqrt(rsum);
    rsum *= 1e6;
    return(rsum);
}

/* Daily avg rms - residuals should be already ~sorted by date */
double resid_rms_daily(int npts, struct residual *r, int weight) {
    int i;
    double rsum=0.0, wsum=0.0, drsum=0.0, dwsum=0.0, w;
    w = weight ? r[0].weight : 1.0; 
    drsum = w*r[0].res_sec;
    dwsum = w;
    for (i=1; i<npts; i++) {
        if (fabs(r[i].toa-r[i-1].toa)>1.0) { 
            drsum /= dwsum;
            rsum += dwsum*drsum*drsum;
            wsum += dwsum;
            drsum=0.0;
            dwsum=0.0;
        } 
        w = weight ? r[i].weight : 1.0; 
        drsum += w*r[i].res_sec;
        dwsum += w;
    }
    drsum /= dwsum;
    rsum += dwsum*drsum*drsum;
    wsum += dwsum;
    rsum /= wsum;
    rsum = sqrt(rsum);
    rsum *= 1e6;
    return(rsum);
}

int main(int argc, char *argv[]) {
    int rv;

    /* Parse opts */
    static struct option long_opts[] = {
        {"help",       0, NULL, 'h'},
        {"stats",      0, NULL, 's'},
        {"bands",      0, NULL, 'b'},
        {"mjd",        0, NULL, 'm'},
        {"res_phase",  0, NULL, 'p'},
        {"res_sec",    0, NULL, 't'},
        {"res_us",     0, NULL, 'r'},
        {"ophase",     0, NULL, 'o'},
        {"rf",         0, NULL, 'f'},
        {"weight",     0, NULL, 'w'},
        {"err",        0, NULL, 'e'},
        {"prefit_sec", 0, NULL, 'i'},
        {"ddm",        0, NULL, 'd'},
        {0,0,0,0}
    };
    int opt, opti;
    char outputs[MAX_OUTS];
    int nout=0;
    outputs[0]='\0';
    int do_stat=0, do_band=0;
    while ((opt=getopt_long(argc,argv,"hsbmptrofweid",long_opts,&opti))!=-1) {
        switch (opt) {
            case 'm':
            case 'p':
            case 't':
            case 'r':
            case 'o':
            case 'f':
            case 'w':
            case 'e':
            case 'i':
            case 'd':
                outputs[nout] = opt;
                nout++;
                if (nout==MAX_OUTS) { 
                    fprintf(stderr, "Too many options.\n");
                    usage();
                    exit(1);
                }
                outputs[nout] = '\0';
                break;
            case 's':
                do_stat=1;
                break;
            case 'b':
                do_band=1;
                break;
            case 'h':
            default:
                usage();
                exit(0);
                break;

        }
    }

    /* Fill default options if none given */
    if (nout==0) { sprintf(outputs, "mfreo"); nout=5; }

    /* Use default filename if none given */
    char fname[256];
    if (optind>=argc) {
        sprintf(fname, "resid2.tmp");
    } else {
        sprintf(fname, "%s", argv[optind]);
    }

    FILE *rf=NULL;
    struct residual *r=NULL;
    int npts;
    rf = fopen(fname, "r");
    if (rf==NULL) {
        fprintf(stderr, "Error opening %s\n", fname);
        if (optind>=argc) { usage(); }
        exit(1);
    }
    rv = read_resid(rf, &r, &npts);
    if (rv<0) { 
        fprintf(stderr, "Error parsing %s\n", fname);
        fclose(rf);
        exit(0);
    }
    fclose(rf);

    /* Print some summary statistics at top */
    int i,j;
    if (do_stat) {

        /* Basic info */
        printf("# File: %s\n", fname);
        printf("# Ntoa: %d\n", npts);

        /* Compute anything else and put it here */
        printf("#  RMS(raw): %.4f\n", resid_rms(npts, r, 0));
        printf("# WRMS(raw): %.4f\n", resid_rms(npts, r, 1));
        printf("#  RMS(day): %.4f\n", resid_rms_daily(npts, r, 0));
        printf("# WRMS(day): %.4f\n", resid_rms_daily(npts, r, 1));

        /* Last line, describes outputs */
        printf("# Data: ");
        for (i=0; i<nout; i++) {
            switch (outputs[i]) {
                case 'm':
                    printf("MJD");
                    break;
                case 'p':
                    printf("Resid(turns)");
                    break;
                case 't':
                    printf("Resid(s)");
                    break;
                case 'r':
                    printf("Resid(us)");
                    break;
                case 'o':
                    printf("Ophase(turns)");
                    break;
                case 'f':
                    printf("RF(MHz)");
                    break;
                case 'w':
                    printf("Weight");
                    break;
                case 'e':
                    printf("Err(us)");
                    break;
                case 'i':
                    printf("Prefit(s)");
                    break;
                case 'd':
                    printf("DDM");
                    break;
            }
            if (i==nout-1) { printf("\n"); } else { printf(" "); }
        }
    }

    /* Print outputs */
    float last_rf=0.0;
    for (i=0; i<npts; i++) {
        if ((fabs(r[i].rf_bary-last_rf)>100.0) && (i>0) && do_band) 
            printf("\n\n"); 
        for (j=0; j<nout; j++) {
            switch (outputs[j]) {
                case 'm':
                    printf("%15.9f", r[i].toa);
                    break;
                case 'p':
                    printf("%+.8e", r[i].res_phase);
                    break;
                case 't':
                    printf("%+.8e", r[i].res_sec);
                    break;
                case 'r':
                    printf("%+.8e", r[i].res_sec*1e6);
                    break;
                case 'o':
                    printf("%.8f", r[i].ophase);
                    break;
                case 'f':
                    printf("%9.4f", r[i].rf_bary);
                    break;
                case 'w':
                    printf("%.4e", r[i].weight);
                    break;
                case 'e':
                    printf("%6.3e", r[i].err_us);
                    break;
                case 'i':
                    printf("%+.8e", r[i].prefit_sec);
                    break;
                case 'd':
                    printf("%9.5f", r[i].ddm);
                    break;
            }
            if (j==nout-1) { printf("\n"); }
            else { printf(" "); }
        }
        last_rf = r[i].rf_bary;
    }

    exit(0);
}

