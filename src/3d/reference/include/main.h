#ifndef __MAIN_HPC__
#define __MAIN_HPC__

#define DEBUG

#ifdef DEBUG
#define DEBUG_PRINT(...) fprintf(stderr, "DEBUG: " __VA_ARGS__);
#else
#define DEBUG_PRINT(...) (void(0));
#endif

#endif
