from src import clean_counts
from src import fastqc
from src import star
from src import subread

if __name__ == '__main__':
    fastqc.run()
    star.run()
    subread.run()
    clean_counts.run()


