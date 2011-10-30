#ifndef IO_H
#define IO_H

#include "state.h"

void write_header(FILE* fp, int n);
void write_frame_data(FILE* fp, int n, particle_t *particles, int* c);

#endif /* IO_H */
