/* @(#)pixfont.h 1.14 88/02/08 SMI */

/*
 * Copyright (c) 1986 by Sun Microsystems, Inc.
 */

#ifndef pixfont_DEFINED
#define pixfont_DEFINED

/*
 * Definition of pixfonts for pixrect library.
 * Include <pixrect/pixrect.h> before this file.
 */

/*
 * A character descriptor contains the pixrect constituting the actual
 * character, the coordinates of the home (top left) of that pixrect
 * relative to the character origin (a point on the baseline near the left
 * of the character), and the distance by which to advance the origin after
 * printing this character.
 */
struct pixchar
{
  Pixrect *pc_pr;    /* pixrect for this char */
  struct pr_pos pc_home;    /* home coords relative to left baseline */
  struct pr_pos pc_adv;    /* distance to next char */
};

/*
 * A font descriptor contains the width of a space (intended to be used
 * in computing backspace and tab distances), the distance between consecutive
 * baselines in the absence of any superscripting, subscripting, or similar
 * vertical monkey business, and an array of 256 character descriptors.
 */
typedef struct pixfont
{
  struct pr_size pf_defaultsize; /* default character size */
  struct pixchar pf_char[256];
}
Pixfont;

/* structured text macros */
#ifndef lint
#define prs_text(prpos, op, pf, str) \
 pr_text((prpos).pr, (prpos).pos.x, (prpos).pos.y, pf, str)

#define prs_ttext(prpos, op, pf, str) \
 pr_ttext((prpos).pr, (prpos).pos.x, (prpos).pos.y, pf, str)
#endif /*lint*/

Pixfont *pf_open();
Pixfont *pf_open_private();
Pixfont *pf_default();
struct pr_size pf_textbatch();
struct pr_size pf_textwidth();

#define PIXFONT Pixfont

#endif /*pixfont_DEFINED*/
