from src import clean_counts, diff_gene_expression, fastqc, get_example_fastq, star, subread

if __name__ == '__main__':
    get_example_fastq.run()
    fastqc.run()
    star.run()
    subread.run()
    clean_counts.run()
    diff_gene_expression.run()


