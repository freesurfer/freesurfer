/**
 * @brief Routine to check .license file
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef Darwin
#include <gnu/libc-version.h>
#endif

#include <const.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#ifndef NO_FIPS_SUPPORT
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#endif

#include "chklc.h"
#include "diag.h"

#define MAX_KEY_LEN 1024

static const char *errmsg =
    "--------------------------------------------------------------------------\n"
    "ERROR: FreeSurfer environment FREESURFER_HOME is not defined.\n"
    "  If you are outside the NMR-Martinos Center, please set this\n"
    "  variable to the location where you installed FreeSurfer.\n"
    "  If you are inside the NMR-Martinos Center, please source\n"
    "  the standard environment. If you need to install FreeSurfer,\n"
    "  go to: http://surfer.nmr.mgh.harvard.edu\n"
    "--------------------------------------------------------------------------\n";

static const char *licmsg =
    "--------------------------------------------------------------------------\n"
    "ERROR: FreeSurfer license file %s not found.\n"
    "  If you are outside the NMR-Martinos Center,\n"
    "  go to http://surfer.nmr.mgh.harvard.edu/registration.html to \n"
    "  get a valid license file (it's free).\n"
    "  If you are inside the NMR-Martinos Center,\n"
    "  make sure to source the standard environment.\n"
    "  A path to an alternative license file can also be\n"
    "  specified with the FS_LICENSE environmental variable.\n"
    "--------------------------------------------------------------------------\n";

static const char *licmsg2 =
    "--------------------------------------------------------------------------\n"
    "ERROR: Invalid FreeSurfer license key found in license file %s\n"
    "  If you are outside the NMR-Martinos Center,\n"
    "  go to http://surfer.nmr.mgh.harvard.edu/registration.html to \n"
    "  get a valid license file (it's free).\n"
    "  If you are inside the NMR-Martinos Center,\n"
    "  make sure to source the standard environment.\n"
    "--------------------------------------------------------------------------\n";

static const char *permission_msg =
    "---------------------------------------------------------------------------\n"
    "ERROR: FreeSurfer license file %s exists "
    "but you do not have read permission.\n"
    "Try running:\n\n"
    "  chmod a+r %s\n"
    "---------------------------------------------------------------------------\n";

static const char *isdir_msg =
    "---------------------------------------------------------------------------\n"
    "ERROR: FS_LICENSE environment variable points to a folder not a file\n"
    "---------------------------------------------------------------------------\n";

#ifndef NO_FIPS_SUPPORT
// 256 bit key
static unsigned char aes_256_cbc_key[] = { 0x3c, 0x3f, 0x78, 0x6d, 0x6c, 0x20, 0x76, 0x65,
			                   0x72, 0x73, 0x69, 0x6f, 0x6e, 0x3d, 0x22, 0x31,
			                   0x2e, 0x30, 0x22, 0x20, 0x65, 0x6e, 0x63, 0x6f,
                                           0x64, 0x69, 0x6e, 0x67, 0x3d, 0x22, 0x49, 0x53
                                         };
// 128 bit IV
static unsigned char aes_256_cbc_iv[]  = { 0x6f, 0x73, 0x69, 0x74, 0x69, 0x6f, 0x6e, 0x61,
			                   0x6c, 0x2a, 0x20, 0x2c, 0x20, 0x72, 0x65, 0x71
                                         };

static int __decrypt_openssl(unsigned char *ciphertext, int ciphertext_len, unsigned char *key,
            unsigned char *iv, unsigned char *plaintext);	    
static int __handleErrors_openssl();
#endif

// freeview will pass a msg buffer to chklc() call
int chklc(char *msg)
{
  static int first_time = 1;
  
  char str[STRLEN];
  if (msg != NULL)
    sprintf(str, "S%sER%sRONT%sOR", "URF", "_F", "DO");
  else
    sprintf(str, "S%sER%sIDE%sOR", "URF", "_S", "DO");
  
  if (getenv(str) != NULL)
    return 1;

  char dirname[STRLEN] = {'\0'};
  char *cp = getenv("FREESURFER_HOME");
  if (cp == NULL) {
    fprintf(stderr, "%s", errmsg);
#ifdef Darwin
    fprintf(stderr, "\n");
    fprintf(stderr, "%s", "Attempting to use the /Applications/freesurfer directory.\n");
    strncpy(dirname, "/Applications/freesurfer", STRLEN);
#else
    if (msg != NULL)
    {
      sprintf(msg, "%s", errmsg);
      return 0;
    }
    exit(-1);
#endif
  }
  else {
    strncpy(dirname, cp, STRLEN-1);
  }

  char lfilename[STRLEN] = {'\0'};
  FILE *lfile = NULL;
  
  // check if alternative license path is provided:
  char *alt = getenv("FS_LICENSE");
  if (alt != NULL) {
    strncpy(lfilename, alt, 511);	// leave a nul on the end
    if (Gdiag_no > 0 && first_time) printf("Trying license file %s\n", lfilename);
    lfile = fopen(lfilename, "r");
    if (lfile == NULL) {
      if (errno == EACCES) {
        fprintf(stderr, permission_msg, lfilename, lfilename);
	if (msg != NULL)
        {
	  sprintf(msg, permission_msg, lfilename, lfilename);
	  return 0;
        }
        exit(-1);
      }
      fprintf(stderr, licmsg, lfilename);
      if (msg != NULL)
      {
	sprintf(msg, licmsg, lfilename);
	return 0;
      }
      exit(-1);
    }
    // make sure that the path is not a directory
    struct stat path_stat;
    stat(alt, &path_stat);
    if S_ISDIR(path_stat.st_mode) {;
      puts(isdir_msg);
      exit(-1);
    }
  }

  // check for license in FREESURFER_HOME:
  if (lfile == NULL) {
    auto cx = snprintf(lfilename, 511, "%s/lic%s", dirname, "ense.txt");
    if( (cx<0) || (cx>511) ) {
      std::cerr << __FUNCTION__
		<< ": snprintf returned error on line "
		<< __LINE__ << std::endl;
    }
    if (Gdiag_no > 0 && first_time) printf("Trying license file %s\n", lfilename);
    lfile = fopen(lfilename, "r");
  }
  if (lfile == NULL) {
    if (errno == EACCES) {
      fprintf(stderr, permission_msg, lfilename, lfilename);
      if (msg != NULL)
      {
        sprintf(msg, permission_msg, lfilename, lfilename);
	return 0;
      }
      exit(-1);
    }
    auto cx = snprintf(lfilename, 511, "%s/.lic%s", dirname, "ense");
    if( (cx<0) || (cx>511) ) {
      std::cerr << __FUNCTION__
		<< ": snprintf returned error on line "
		<< __LINE__ << std::endl;
    }
    if (Gdiag_no > 0 && first_time) printf("Now trying license file %s\n", lfilename);
    lfile = fopen(lfilename, "r");
  }
  if (lfile == NULL) {
    if (errno == EACCES) {
      fprintf(stderr, permission_msg, lfilename, lfilename);
      if (msg != NULL)
      {
        sprintf(msg, permission_msg, lfilename, lfilename);
	return 0;
      }      
      exit(-1);
    }
    fprintf(stderr, licmsg, lfilename);
    if (msg != NULL)
    {
      sprintf(msg, licmsg, lfilename);
      return 0;
    }
    exit(-1);
  }

  char email[MAX_KEY_LEN/2] = {'\0'}, magic[MAX_KEY_LEN/2] = {'\0'};
  char key[MAX_KEY_LEN] = {'\0'}, key2[MAX_KEY_LEN] = {'\0'}, key3[MAX_KEY_LEN] = {'\0'};
  if (fscanf(lfile, "%s\n%s\n%s\n%s\n%s\n", email, magic, key, key2, key3) < 4) {
    fprintf(stderr, "error parsing license file, at least 4 values expected\n");
    fclose(lfile);
    fprintf(stderr, licmsg2, lfilename);
    if (msg != NULL)
    {
      sprintf(msg, licmsg2, lfilename);
      return 0;
    }
    exit(-1);
  }
  fclose(lfile);
  
  char gkey[MAX_KEY_LEN] = {'\0'};  
  sprintf(gkey, "%s.%s", email, magic);

  if (Gdiag_no > 0 && first_time) {
    printf("email %s\n", email);
    printf("magic %s\n", magic);
    printf("key   %s\n", key);
    printf("key2  %s\n", key2);
    printf("key3  %s\n", key3);
    printf("gkey  %s\n", gkey);
  }

  // This code is meant to provide backwards compatibility
  // of freesurfer license checking. Unfortunately previous
  // freesurfer license keys used an improper salt value of
  // "*C" which is not valid because it should be purely
  // alpha-numeric. New license files have a 4th line with
  // a key generated using a proper salt.
  char *crypt_gkey = NULL;
#ifndef NO_FIPS_SUPPORT  
  if (strcmp(key3, "") != 0) { // 5 line license file
    if (Gdiag_no > 0)
      printf("[DEBUG] chklc() 5 line license file %s (len=%lu, key3=%s)\n", lfilename, strlen(key3), key3);

    int decryptedkey_len = 0;
    
    // base64 decode the key
    unsigned char decoded_encryptedkey[MAX_KEY_LEN] = {'\0'};
    int decoded_encryptedkey_len = EVP_DecodeBlock(decoded_encryptedkey, (unsigned char*)key3, strlen(key3));
    if (decoded_encryptedkey_len < 0)
    {
      if (Gdiag_no > 0)
        printf("[DEBUG] key3 len = %lu, decoded_encryptedkey_len = %d\n", strlen(key3), decoded_encryptedkey_len);

      printf("ERROR: EVP_DecodeBlock() failed with 5-line file (%s)\n", lfilename);
      fprintf(stderr, licmsg2, lfilename);
      if (msg != NULL)
      {
        sprintf(msg, licmsg2, lfilename);
        return 0;
      }
      exit(-1);
    }

    // Decrypt the decoded_encryptedkey
    unsigned char decryptedkey[MAX_KEY_LEN] = {'\0'};
    if (Gdiag_no > 0)
      printf("[DEBUG] key3 len = %lu, decoded_encryptedkey_len = %d, strlen((char*)decoded_encryptedkey) = %lu\n",
	     strlen(key3), decoded_encryptedkey_len, strlen((char*)decoded_encryptedkey));
    
    // pass decoded_encryptedkey_len to the call instead of strlen((char*)decoded_encryptedkey)
    //int decryptedkey_len = __decrypt_openssl(decoded_encryptedkey, strlen((char*)decoded_encryptedkey), aes_256_cbc_key, aes_256_cbc_iv, decryptedkey);
    decryptedkey_len = __decrypt_openssl(decoded_encryptedkey, decoded_encryptedkey_len, aes_256_cbc_key, aes_256_cbc_iv, decryptedkey);
    if (decryptedkey_len < 0)
    {
      printf("ERROR: __decrypt_openssl() failed with 5-line file (%s)\n", lfilename);
      fprintf(stderr, licmsg2, lfilename);
      if (msg != NULL)
      {
        sprintf(msg, licmsg2, lfilename);
        return 0;
      }
      exit(1);
    }
    
    decryptedkey[decryptedkey_len] = '\0';
    if (Gdiag_no > 0)
      printf("[DEBUG] chklc() decrypted key: %s (%d) \n", decryptedkey, decryptedkey_len);

    memset(key, 0, MAX_KEY_LEN);
    memcpy(key, gkey, strlen(gkey));
    
    crypt_gkey = (char*)malloc(decryptedkey_len+1);
    memset(crypt_gkey, 0, decryptedkey_len+1);
    memcpy(crypt_gkey, decryptedkey, decryptedkey_len);
    crypt_gkey[strlen(key)] = '\0';

    if (Gdiag_no > 0)
      printf("[DEBUG] key = <%s> (%lu), crypt_gkey = <%s> (%d, strlen=%lu)\n", key, strlen(key), crypt_gkey, decryptedkey_len, strlen(crypt_gkey));
  }
  else if (strcmp(key2, "") != 0) {
#else  
  if (strcmp(key2, "") != 0) {
#endif    
    // We have a 4 line license file.
    if (Gdiag_no > 0)
      printf("[DEBUG] chklc() 4 line license file %s\n", lfilename);
    if (Gdiag_no > 0 && first_time) printf("4 line license file\n");
    strcpy(key, key2);
    crypt_gkey = crypt(gkey, "FS");
    if (crypt_gkey == NULL) {
      printf("ERROR: crypt() returned null with 4-line file (%s)\n", lfilename);
      if (msg != NULL)
      {
	sprintf(msg, "ERROR: crypt() returned null with 4-line file (%s)\n", lfilename);
	return 0;
      }
      exit(1);
    }
  }
  else {
    // We have a 3 line license file.
    if (Gdiag_no > 0 && first_time) printf("3 line license file\n");
#ifdef Darwin
    // On Darwin systems the key produced with a salt of '*C'
    // is different than that produced on Linux. So to be backwards
    // compatible we must avoid the call to crypt.
    crypt_gkey = key;
#else
    cmp_glib_version();
    crypt_gkey = crypt(gkey, "*C");
#endif
  }

  if (Gdiag_no > 0 && first_time) printf("crypt_gkey %s\n", crypt_gkey);

  if (memcmp(key, crypt_gkey, strlen(key)) != 0) {
    fprintf(stderr, licmsg2, lfilename);
    if (msg != NULL)
    {
      sprintf(msg, licmsg2, lfilename);
      return 0;
    }
    exit(-1);
  }

#ifndef NO_FIPS_SUPPORT  
  if (strcmp(key3, "") != 0)
    free(crypt_gkey);
#endif  
  
  if (Gdiag_no > 0 && first_time) printf("chklc() done\n");
  first_time = 0;

  return 1;
}

#ifndef Darwin
void cmp_glib_version(void)
{
  int i;
  const char *GNU_LIBC_VERSION_MAX = "2.15";
  int glibc_max[2], glibc_current[2];
  static const char *new_license_msg =
      "--------------------------------------------------------------------------\n"
      "GNU libc version: %s\n"
      "ERROR: Systems running GNU glibc version greater than %s\n"
      "  require a newly formatted license file (it's free). Please\n"
      "  download a new one from the following page:\n"
      "  http://surfer.nmr.mgh.harvard.edu/registration.html\n"
      "--------------------------------------------------------------------------\n";

  sscanf(GNU_LIBC_VERSION_MAX, "%d.%d", &glibc_max[0], &glibc_max[1]);
  sscanf(gnu_get_libc_version(), "%d.%d", &glibc_current[0], &glibc_current[1]);

  for (i = 0; i < 2; i++) {
    if (glibc_current[i] > glibc_max[i]) {
      printf(new_license_msg, gnu_get_libc_version(), GNU_LIBC_VERSION_MAX);
      exit(-1);
    }
  }
}
#endif

#ifndef NO_FIPS_SUPPORT
// return -1 for errors
static int __handleErrors_openssl()
{
  ERR_print_errors_fp(stderr);
  //abort();
  return -1;
}

// return the decrypted text length
// return -1 for errors
static int __decrypt_openssl(unsigned char *ciphertext, int ciphertext_len, unsigned char *key,
                             unsigned char *iv, unsigned char *plaintext)
{
  EVP_CIPHER_CTX *ctx;

  /* Create and initialise the context */
  if (!(ctx = EVP_CIPHER_CTX_new()))
    return __handleErrors_openssl();

  /*
   * Initialise the decryption operation
   */
  if (1 != EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv))
    return __handleErrors_openssl();

  // get cipher block size, it is 128 bits (16 bytes) for NIST AES standard
  int cipher_block_size = EVP_CIPHER_CTX_block_size(ctx);
  if (Gdiag_no > 0)
    printf("[DEBUG] cipher_block_size= %d\n", cipher_block_size);

  // disable padding
  EVP_CIPHER_CTX_set_padding(ctx, 0);

  // make the text to be decrypted multiple of cipher block size
  int todecrypttext_len = (ciphertext_len%cipher_block_size == 0) ? ciphertext_len : (ciphertext_len/cipher_block_size + 1)*cipher_block_size;
  unsigned char todecrypttext[todecrypttext_len];
  memset(todecrypttext, 0, todecrypttext_len);
  memcpy(todecrypttext, ciphertext, ciphertext_len);
  if (Gdiag_no > 0)
    printf("[DEBUG] ciphertext_len = %d, todecrypttext_len = %d\n", ciphertext_len, todecrypttext_len);
    
  /*
   * Provide the message to be decrypted, and obtain the plaintext output.
   * EVP_DecryptUpdate can be called multiple times if necessary.
   */
  int len;
  if (1 != EVP_DecryptUpdate(ctx, plaintext, &len, todecrypttext, todecrypttext_len))
    return __handleErrors_openssl();

  int plaintext_len = len;

  /*
   * Finalise the decryption.
   * Further plaintext bytes may be written at this stage.
   */
  if (1 != EVP_DecryptFinal_ex(ctx, plaintext + len, &len))
    return __handleErrors_openssl();

  plaintext_len += len;

  /* Clean up */
  EVP_CIPHER_CTX_free(ctx);

  return plaintext_len;
}
#endif
