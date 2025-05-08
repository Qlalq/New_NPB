static void print_timers(void) {
    int i;
    char *tstrings[] = { "          total ",
			 "          setup ", 
			 "            fft ", 
			 "         evolve ", 
			 "       checksum ", 
			 "         fftlow ", 
			 "        fftcopy " };
    for (i = 0; i < T_MAX; i++) {
	if (timer_read(i) != 0.0) {
            printf("timer %2d(%16s( :%10.6f\n", i, tstrings[i], timer_read(i));
	}
    }
}