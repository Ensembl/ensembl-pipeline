/* Silly little program to rapidly extract subsequences from single
   entry FASTA files.

   Note that it makes the following assumptions:

   1)  That the FASTA file contains only one entry
   2)  That the sequence lines (except for the last) are all the same length
   3)  That the sequence line length is less than 1023

   The code has been tested on Tru64 and on Linux/x86, but no other
   platforms.

   Tim Cutts <tjrc@sanger.ac.uk> */

#ifdef linux
#  define _LARGEFILE_SOURCE
#  define _FILE_OFFSET_BITS 64
#endif

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
  char buf[1024];
  FILE *fa;
  off_t start, stop;
  off_t from, to;
  off_t beginning;
  int linewidth;
  char c;

  if (sizeof(off_t) < 8) {
    fprintf(stderr, "WARNING: Not compiled with large file support\n");
  }

  if (argc != 4) {
    fprintf(stderr, "Usage: %s fasta_file start stop\n", argv[0]);
    return 1;
  }

  if ((fa = fopen(argv[1], "r")) == NULL) {
    perror(argv[1]);
    return 2;
  }

  start = atol(argv[2]);
  stop = atol(argv[3]);

  if (start > stop) {
    from = start;
    start = stop;
    stop = from;
  }

  if (fgets(buf, 1023, fa) == NULL) {
    perror(argv[1]);
    return 2;
  }

  if (buf[0] != '>') {
    fprintf(stderr, "%s: Not a fasta file\n", argv[1]);
    return 1;
  }

  beginning = ftello(fa);

  if (fgets(buf, 1023, fa) == NULL) {
    perror(argv[1]);
    return 2;
  }

  linewidth = strlen(buf) - 1;

  if (buf[linewidth] != '\n') {
    fprintf(stderr, "ERROR: FASTA sequence lines are too wide\n");
    return 1;
  }

  from = beginning + start + ((start - 1)/ linewidth) - 1;
  to = beginning + stop + ((stop - 1)/ linewidth) - 1;

  if (fseeko(fa, to, SEEK_SET) < 0) {
    perror("Seek error");
    return 1;
  }

  fseeko(fa, from, SEEK_SET);

  while (from++ <= to) {
    c = fgetc(fa);
    if (c != '\n') {
      if (c == '>' || c == EOF) {
	fprintf(stderr,
		"WARNING: %d is past the end of the sequence\n", stop);
	break;
      }
      putc(c, stdout);
    }
  }

  fclose(fa);

  putc('\n', stdout);

  return 0;
}
