CC = cc
CFLAGS = -Wall -Wextra -O2 -std=c11

.PHONY: benchmark clean

benchmark: benchmark_dh

benchmark_dh: benchmark_dh.c
	$(CC) $(CFLAGS) benchmark_dh.c -o benchmark_dh

clean:
	rm -f benchmark_dh
