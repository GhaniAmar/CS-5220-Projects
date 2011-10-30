#ifndef IO_H
#define IO_H

#include "state.h"

void write_header(FILE* fp, int n);
void write_frame_data(FILE* fp, sim_state_t *state, int *c);

#endif /* IO_H */
