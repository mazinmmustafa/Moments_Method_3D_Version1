
FiguresGP = $(wildcard *.gp)
FiguresPDF = $(patsubst %.gp, %.pdf, $(FiguresGP))

all: clean $(FiguresPDF)

clean:
	$(RM) -r $(FiguresPDF)

%.pdf: %.gp
	gnuplot $^
	pdflatex $(patsubst %.pdf, %.tex, $@) > dump.txt