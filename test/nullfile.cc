#include <cstdio>

int main(int argc, char ** argv)
{

    FILE * cdfs_output_fh = fopen(argv[1], "w");
    if (cdfs_output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open cdfs_output_file %s\n",
                argv[1]);
    }
    return 0;
}
