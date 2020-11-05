GFF = HT.Gene.gff
ZIP = HT.KYGene.gff3.zip
URL = http://ghost.zool.kyoto-u.ac.jp/datas/$(ZIP)

install: 
	Rscript -e "devtools::document()"
	R CMD INSTALL .

data/atacCiona.db: data data-raw/$(GFF)
	Rscript data-raw/writeDB.R
	Rscript data-raw/DA.R

data-raw/$(GFF): data
	wget -U firefox $(URL)
	unzip -o $(ZIP)
	mv $(GFF) data-raw

data:
	mkdir -p data

clean:
	rm -f $(ZIP)
