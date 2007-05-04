/*$Header: /space/repo/1/dev/dev/talairach_avi/ANALYZE.h,v 1.1 2007/05/04 22:33:58 nicks Exp $*/
/*$Log: ANALYZE.h,v $
/*Revision 1.1  2007/05/04 22:33:58  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.1  2005/12/15  22:54:50  avi
 * Initial revision
 **/
 
/***************************************/
/* struct dsr--the ANALYZE .hdr struct */
/***************************************/
struct header_key{
        int sizeof_hdr;                 /*required--byte size of header file*/
        char data_type[10];
        char db_name[18];
        int extents;                    /*required--16384*/
        short int session_error;
        char regular;                   /*required--'r'=regular*/
        char hkey_un0;
};
struct image_dimension{
        short int dim[8];               /*required*/       
        short int unused8;
        short int unused9;
        short int unused10;
        short int unused11;
        short int unused12;
        short int unused13;
        short int unused14;
        short int datatype;             /*required: 0=unk,1=1 bit/pix,2=8bits,4=16 bits*/
                                        /*8=32 bits (signed int),16=32 bits (floating pt)*/
                                        /*32=64 bits (2 floats),64=64 bits (double)     */
        short int bitpix;               /*bits/pixel*/
        short int dim_un0;
        float pixdim[8];                /*real world values of dimensions mm ms*/
        float funused8;
        float funused9;
        float funused10;
        float funused11;
        float funused12;
        float funused13;
        float compressed;
        float verified;
        int glmax,glmin;                /*required*/
};
struct data_history{
        char descrip[80];               /*Will be displayed when loading*/
        char aux_file[24];
        char orient;
        char originator[10];
        char generated[10];
        char scannum[10];
        char patient_id[10];
        char exp_date[10];
        char exp_time[10];
        char hist_un0[3];
        int views;
        int vols_added;
        int start_field;
        int field_skip;
        int omax,omin;
        int smax,smin;
};
struct dsr{
        struct header_key hk;
        struct image_dimension dime;
        struct data_history hist;
};
