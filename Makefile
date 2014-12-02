CFLAGS = -O3 -Wall -Wextra -fno-strict-aliasing

all: daligner daligner_p HPCdaligner \
     LAsort LAmerge LAsplit LAcat LAshow LAcheck LA4Falcon DB2Falcon DB.so \

daligner: daligner.c filter.c filter.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner daligner.c filter.c align.c DB.c QV.c -lpthread -lm

daligner_p: daligner.c filter_p.c filter.h align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o daligner_p daligner.c filter_p.c align.c DB.c QV.c -lpthread -lm

HPCdaligner: HPCdaligner.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o HPCdaligner HPCdaligner.c DB.c QV.c -lm

LAsort: LAsort.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsort LAsort.c DB.c QV.c -lm

LAmerge: LAmerge.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAmerge LAmerge.c DB.c QV.c -lm

LAshow: LAshow.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAshow LAshow.c align.c DB.c QV.c -lm

LA4Falcon: LA4Falcon.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LA4Falcon LA4Falcon.c align.c DB.c QV.c -lm

LAcat: LAcat.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcat LAcat.c DB.c QV.c -lm

LAsplit: LAsplit.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAsplit LAsplit.c DB.c QV.c -lm

LAcheck: LAcheck.c align.c align.h DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o LAcheck LAcheck.c align.c DB.c QV.c -lm

DB2Falcon: DB2Falcon.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DB2Falcon DB2Falcon.c DB.c QV.c -lm

DB.so: DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -shared -fPIC -o DB.so DB.c DB.h QV.c -lm


clean:
	rm -f daligner HPCdaligner
	rm -f LAsort LAmerge LAshow LAsplit LAcat LAcheck DB.so LA4Falcon
	rm -f DB2Falcon

install:
	cp daligner HPCdaligner ~/bin
	cp LAsort LAmerge LAshow LAsplit LAcat LAcheck ~/bin

package:
	make clean
	tar -zcf daligner.tar.gz README *.h *.c Makefile

