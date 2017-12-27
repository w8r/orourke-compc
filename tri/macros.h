char *malloc();

#define NEW(p, type) if ((p=(type *) malloc (sizeof(type))) == NULL) {\
      printf ("NEW: Out of Memory!\n");\
      exit(EXIT_FAILURE);\
    }

#define ADD( head, p )  if ( head )  { \
       p->next = head; \
       p->prev = head->prev; \
       head->prev = p; \
       p->prev->next = p; \
    } \
    else { \
       head = p; \
       head->next = head->prev = p; \
    }
