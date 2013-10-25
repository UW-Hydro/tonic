#!/usr/local/bin/python

def replace(template_path, outfile_path, line_pattern, replace_str):
    #Create temp file
    out_file = open(outfile_path,'w')
    template_file = open(template_path)
    for line in template_file:
        out_file.write(line.replace(line_pattern, replace_str))
    #close files
    template_file.close()
    out_file.close()

if __name__ == "__main__":
    main()


