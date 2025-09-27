cat \
    <(grep active_hor hg002v1.0.1.cenSatv2.0.noheader.bed | \
    bedtools merge -i - -d 5000000 | \
    # https://www.gnu.org/software/gawk/manual/html_node/Round-Function.html
    awk -v OFS="\t" -v OFMT='%f' '
    function round(x,   ival, aval, fraction) {
        ival = int(x)    # integer part, int() truncates

        # see if fractional part
        if (ival == x)   # no fraction
            return ival   # ensure no decimals

        if (x < 0) {
            aval = -x     # absolute value
            ival = int(aval)
            fraction = aval - ival
            if (fraction >= .5)
                return int(x) - 1   # -2.5 --> -3
            else
                return int(x)       # -2.3 --> -2
        } else {
            fraction = x - ival
            if (fraction >= .5)
                return ival + 1
            else
                return ival
        }
    }
    {
        midpt=round((($3-$2)/2) + $2);
        print $1, $2, midpt, ".", "acen";
        print $1, midpt, $3, ".", "acen"
    }') \
    cytoBand.hg002v1.0.bed | \
    sort -k1,1 -k2,2n > /tmp/cytoBand.hg002v1.0.sorted.bed

bedtools subtract \
    -a <(grep -v <(printf "chrEBV\nchrM\n") results/curated/data/asm/v1.0.1.fa.gz.fai | awk -v OFS="\t" '{ print $1, 0, $2}') \
    -b /tmp/cytoBand.hg002v1.0.sorted.bed | awk -v OFS="\t" '{print $0, ".", "none"}' >> /tmp/cytoBand.hg002v1.0.sorted.bed

sort -k1,1 -k2,2n /tmp/cytoBand.hg002v1.0.sorted.bed > results/cytoBand.hg002v1.0.sorted.bed
rm /tmp/cytoBand.hg002v1.0.sorted.bed
