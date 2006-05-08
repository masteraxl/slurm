##*****************************************************************************
#  $Id$
##*****************************************************************************
#  AUTHOR:
#    Morris Jette <jette1@llnl.gov>
#
#  SYNOPSIS:
#    X_AC_PAM
#
#  DESCRIPTION:
#    Test for PAM (Pluggable Authentication Module) support.
#
#  WARNINGS:
#    This macro must be placed after AC_PROG_CC or equivalent.
##*****************************************************************************

AC_DEFUN([X_AC_PAM], [
  AC_MSG_CHECKING([whether to enable PAM support])
  AC_ARG_ENABLE(
    [pam],
    AC_HELP_STRING([--enable-pam], [enable PAM (Pluggable Authentication Modules) support]),
    [ case "$enableval" in
        yes) x_ac_pam=yes ;;
         no) x_ac_pam=no ;;
          *) AC_MSG_RESULT([doh!])
             AC_MSG_ERROR([bad value "$enableval" for --enable-pam]) ;;
      esac
    ],
    [x_ac_pam=yes]
  )

  if test "$x_ac_pam" = yes; then
    AC_MSG_RESULT([yes])
    AC_CHECK_LIB([pam],
        [pam_get_user],
        [ac_have_pam=yes; PAM_LIBS="-lpam"])

    AC_CHECK_LIB([pam_misc],
        [misc_conv],
        [ac_have_pam_misc=yes; PAM_LIBS="$PAM_LIBS -lpam_misc"])

    AC_SUBST(PAM_LIBS)
    if test "x$ac_have_pam" = "xyes" -a "x$ac_have_pam_misc" = "xyes"; then
      AC_DEFINE(HAVE_PAM,, [define if you have the PAM library])
    else
      AC_MSG_WARN([Unable to locate PAM libraries])
    fi
  else
    AC_MSG_RESULT([no])
  fi

  AM_CONDITIONAL(HAVE_PAM,
      test "x$x_ac_pam" = "xyes" -a "x$ac_have_pam" = "xyes" -a "x$ac_have_pam_misc" = "xyes")

])

