import tabula 
install tabula
file = "data/Pam_chip/STK-144-PamChip-87102.pdf"
output = tabula.convert_into(file, "data/Pam_chip/Pam_converted.csv", output_format="csv", lattice=True, stream=False,  pages="all" )
