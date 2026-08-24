dir = getDirectory( "Choose the Directory" );
list = getFileList( dir );
for ( i=0; i<list.length; i++ ) {
    open( dir + list[i] );
}
run("Concatenate...", "all_open open");

saveAs("tiff");
