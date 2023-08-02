import pyliftover
# Load the liftover chain file for GRCh38 to GRCh37
lo = pyliftover.LiftOver('hg38', 'hg19')
lo_37_to_38 = pyliftover.LiftOver('hg19', 'hg38')
grch37 = lo_37_to_38.convert_coordinate(
    "chr1",
    227468150
)
if grch37:
    # liftover performed successfully
    (
        chromosome_grch37,
        start_grch37,
        _,
        _
    ) = grch37[0]
    print(start_grch37)