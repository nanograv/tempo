/*  $Id$   */

#include <stdio.h>

/*  FORTRAN:  fd = close(filedes)      */
close_(filedes)
int *filedes;
{
return(close(*filedes));
}
/*  FORTRAN:  fd = open(filnam,mode)  */
open_(filnam,mode)
char filnam[];
int *mode;
{
  return(open(filnam,*mode));
}
/* FORTRAN:  fd = creat(filnam,mode) */
creat_(filnam,mode)
char filnam[];
int *mode;
{
  return(creat(filnam,*mode));
}
/* FORTRAN:  nread = read(fd,buf,n) */
read_(fd,buf,n)
int *fd,*n;
char buf[];
{
  return(read(*fd,buf,*n));
}
/* FORTRAN:  nwrt = write(fd,buf,n) */
write_(fd,buf,n)
int *fd,*n;
char buf[];
{
  return(write(*fd,buf,*n));
}
/* FORTRAN: ns = lseek(fd,offset,origin) */
lseek_(fd,offset,origin)
int *fd,*offset,*origin;
{
  return(lseek(*fd,*offset,*origin));
}
/* times(2) */
times_(buf)
int buf[];
{
  return (times(buf));
}
/* returns 32 bit random numbers to FORTRAN program */
iran_(xsubi)
unsigned short xsubi[];
{
  return (jrand48(xsubi));
}
exit_(n)
int *n;
{
  printf("\n\n");
  exit(*n);
}

cprnt_(string,n)
char string[];
int *n;
{
  int i;

  for(i=0;i< *n;i++)
    printf("%s",string);
  fflush(stdout);
  return;
}

byterev_(i,n)
int *i,*n;
{
int j,k;

for(j=0;j<*n;j++) {
  k  = (*i<<24) | ((*i>>24)&255);
  k  = k | ((*i&65280)<<8);
  *i = k | ((*i>>8)&65280);
  i++;
}
}

dbyterev_(i,n)
int *i,*n;
{
int j,k,l;

for(j=0;j<*n;j++) {
  k = (*i<<24) | ((*i>>24)&255);
  k = k | ((*i&65280)<<8);
  k = k | ((*i>>8)&65280);
  i++;
  l = (*i<<24) | ((*i>>24)&255);
  l = l | ((*i&65280)<<8);
  l = l | ((*i>>8)&65280);
  *(i-1) = l;
  *i = k;
  i++;
}
}

wswap_(i,n)
int *i,*n;
{
int j,k;

for(j=0;j<*n;j++) {
  k = *i;
  *i = *(i+1);
  *(i+1) = k;
  i+=2;
}
}

double hrtime_()
{
  int tv[2],tz[2];
  gettimeofday(tv,tz);
  return(tv[0]+1.e-6*tv[1]);
}

int time_()
{
  return(time(0));
}

int lstat_(char *filename, int *buf)
{
  return(lstat(filename,buf));
}

void usleep_(long int *n)
{
  usleep(*n);
}

/* mallocx and freex are variants on code suggested
   by Richard Dodgson (U.Tasmania)  

   sample use:  if a real*8 array of size N is needed:
      integer ipointer, aoff
      real*8 a(1)
      ipointer = mallocx(a(1),n,8,aoff)

   ... we will compare its use to array b(1:n), which is allocated
       in the normal way:
      real*8 b(n)

   ... now use a(1+aoff) to a(n+aoff) as array references, i.e.,
       add aoff to the index of all references to a(); for example,
       to set the array to 1..n, use:

      do i = 1, n
         a(i+aoff) = i
         b(i) = i
      end do

   ...   when sending the array to a subroutine, use a(1+aoff) in
         the subroutine call [normally, one would write a(1) or a]
      call somesub(a(1+aoff))
      call somesub(b(1))   ...or...
      call somesub(b)

      call freex(ipointer)
*/

void *mallocx_(void *ref, int *nelem, int *size, int *indeks) {
        long nbytes;    /*  long integer necessary on 64-bit systems */
        void *where;
        nbytes = (*nelem)  * (*size) ;
        where = (char *) malloc(nbytes);
        if (where == (char *) NULL)
          { fprintf(stderr,"Out of memory");exit(1); }
        nbytes = (long) ((char *)where - (char *)ref);
        *indeks = nbytes / (*size);
        return where;
}


void freex_(char *(*where)) {
        free(*where);
}



