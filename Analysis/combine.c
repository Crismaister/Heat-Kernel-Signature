#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#define OUTPUT CALCULATE
#define DISTANCE d_euclidian

#define SIMPLE 2
#define SQL 1
#define CALCULATE 3

void readhks (const char *fname, float *hks) {
  FILE *f = fopen(fname, "r");
  char b[100];
  int i;
  assert(f);
  for (i = 0; i < 100; i++) {
    assert(fscanf(f, "%g", hks+i) == 1);
  }
  fclose(f);
}

float d_euclidian (float *a, float *b) {
  float s = 0.0f;
  int i;

  for (i = 0; i < 100; i++) {
    float d = a[i] - b[i];
    s += d*d;
  }
  return sqrtf(s);
}

float d_euclidian_norm_0 (float *a, float *b) {
  float s = 0.0f;
  int i;

  for (i = 0; i < 100; i++) {
    float d = a[i]/a[0] - b[i]/b[0];
    s += d*d;
  }
  return sqrtf(s);
}

float d_euclidian_norm_avg (float *a, float *b) {
  float s = 0.0f;
  int i;

  float avg_a = 0.0f;
  for (i = 0; i < 100; i++)
    avg_a += a[i];
  avg_a /= 100.0;

  float avg_b = 0.0f;
  for (i = 0; i < 100; i++)
    avg_b += b[i];
  avg_b /= 100.0;

  for (i = 0; i < 100; i++) {
    float d = a[i]/avg_a - b[i]/avg_b;
    s += d*d;
  }
  return sqrtf(s);
}

int main (int argc, char **argv) {
  if (argc != 684 + 1) {
    printf("wrong args; usage: ./matrix *.off.signature[.<commit>]\n");
  }
  float hkses[684][100];

  int i, j;

  for (i = 0; i < 684 - 1; i++) {
    assert(strcmp(argv[1 + i], argv[1 + i+1]) < 0);
  }

  for (i = 1; i < argc; i++) {
    readhks(argv[i], hkses[i-1]);
  }

#if OUTPUT == SIMPLE
  for (i = 0; i < 684; i++) {
    for (j = 0; j < 684; j++) {
      printf("%g, ", DISTANCE(hkses[i], hkses[j]));
    }
    printf("\n");
  }
#elif OUTPUT == SQL

  printf("BEGIN TRANSACTION;\n");

  printf("%s",
"CREATE TABLE models ("
"\n  num      int PRIMARY KEY,"
"\n  file     text,"
"\n  model    text,"
"\n  transf   text,"
"\n  strength int);"
"\n"
         );
  printf("%s",
"CREATE TABLE distance ("
"\n  a        int,"
"\n  b        int,"
"\n  d_0      float,"
"\n  d_avg    float,"
"\n  d        float);"
"\n"
         );


  for (i = 0; i < 684; i++) {
    int index = i;
    const char *file = argv[i+1];
    char *fdup;
    char *model;
    char *transf;
    char *strength;
    char *c;
    if (strrchr(file, '/'))
      file = strrchr(file, '/');
    fdup = strdup(file);
    model = fdup;

    c = strchr(model, '.');
    assert(c);
    *c = '\0';
    transf = c+1;

    c = strchr(transf, '.');
    assert(c);
    *c = '\0';
    strength = c+1;

    c = strchr(strength, '.');
    assert(c);
    *c = '\0';

    printf("INSERT INTO models (num, file, model, transf, strength) "
           "VALUES (%d,\"%s\",\"%s\",\"%s\",%s);\n",
           index,
           file,
           model,
           transf,
           strength);
    free(fdup);
  }

  for (i = 0; i < 684; i++) {
    for (j = 0; j < i; j++) {
      printf("INSERT INTO distance(a,b,d,d_0,d_avg) VALUES(%d,%d,%g,%g,%g);\n",
             i, j,
             d_euclidian(hkses[i], hkses[j]),
             d_euclidian_norm_0(hkses[i], hkses[j]),
             d_euclidian_norm_avg(hkses[i], hkses[j])
             );
      printf("INSERT INTO distance(a,b,d,d_0,d_avg) VALUES(%d,%d,%g,%g,%g);\n",
             j, i,
             d_euclidian(hkses[i], hkses[j]),
             d_euclidian_norm_0(hkses[i], hkses[j]),
             d_euclidian_norm_avg(hkses[i], hkses[j])
             );
    }
    printf("INSERT INTO distance(a,b,d,d_0,d_avg) VALUES(%d,%d,%g,%g,%g);\n",
           i, i,
           d_euclidian(hkses[i], hkses[j]),
           d_euclidian_norm_0(hkses[i], hkses[j]),
           d_euclidian_norm_avg(hkses[i], hkses[j])
           );
  }

  printf("CREATE TABLE reference_distances ( d_ref float );");
  float d_ref;
  for (d_ref = 1e-7; d_ref < 2000; d_ref *= 1.05)
    printf("INSERT INTO reference_distances VALUES (%g);\n", d_ref);

  printf("COMMIT TRANSACTION;");
  printf("CREATE VIEW joined AS SELECT"
         " A.model as modelA,"
         " B.model as modelB,"
         " A.transf as transfA,"
         " B.transf as transfB,"
         " A.strength as strengthA,"
         " B.strength as strengthB,"
         " d_avg,"
         " d_0,"
         " d"
         " FROM distance JOIN models A ON a = A.num JOIN models B ON b = B.num;");
#elif OUTPUT == CALCULATE

  float d[684][684];
  int models[684];

  for (i = 0; i < 684; i++) {
    int index = i;
    const char *file = argv[i+1];
    char *fdup;
    char *model;
    char *transf;
    char *strength;
    char *c;
    if (strrchr(file, '/'))
      file = strrchr(file, '/');
    fdup = strdup(file);
    model = fdup;
    assert(sscanf(model, "%i", models+i) == 1);

    /*c = strchr(model, '.');
    assert(c);
    *c = '\0';
    transf = c+1;

    c = strchr(transf, '.');
    assert(c);
    *c = '\0';
    strength = c+1;

    c = strchr(strength, '.');
    assert(c);
    *c = '\0';
    free(fdup);*/
  }

  // good_match[i] == |{p,q | d(p,q) < log(p * 1e6)/log(1.05)}|
  int good_match[502];
  int bad_match[502];

  int relevant_matches = 0;

  for (i = 0; i < 502; i++)
    good_match[i] = bad_match[i] = 0;

  for (i = 0; i < 684; i++) {
    for (j = 0; j < 684; j++) {
      double l;
      int k;
      //printf("%g, ", DISTANCE(hkses[i], hkses[j]));
      d[i][j] = DISTANCE(hkses[i], hkses[j]);
      l = log(d[i][j] * 1e6)/log(1.05);
      if (l < 0) l = 0;
      if (models[i] == models[j])
        for (k = l; k < 500; k++)
          good_match[k]++;
      else
        for (k = l; k < 500; k++)
          bad_match[k]++;
      if (models[i] == models[j])
        relevant_matches++;
    }
  }

  for (i = 0; i < 500; i++) {
    printf("%g %i %i\n", 1e-6 * pow(1.05, i), good_match[i], bad_match[i]);
  }

  // GNUplot:
  // set xlabel "log(Threshold)"; set ylabel "Precision/Recall"; plot "goodbad" using (log($1)):($2/45486) with line title "Recall", "goodbad" using (log($1)):($2/($2+$3)) with line title "Precision"

  // set xlabel "Recall"; set ylabel "Precision"; plot "goodbad" using ($2/45486):($2/($2+$3)) with line title "Precision"
  fprintf(stderr, "relevant: %i\n", relevant_matches);

#else
  #error choose output mode
#endif
}
