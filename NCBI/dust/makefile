CC = cc -O 

.c.o :
	$(CC) -c $<

dust: tust.o dust.o getfa.o
	$(CC) -o dust tust.o dust.o getfa.o
