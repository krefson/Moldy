dnl     Adds to the include-path
dnl
dnl     Autoconf 1.11 should have provided a way to add include path options to
dnl     the cpp-tests.
dnl
AC_DEFUN([CF_INCLUDE_PATH],
[
for cf_path in $1
do
        cf_result="no"
        AC_MSG_CHECKING(for directory $cf_path)
        if test -d $cf_path
        then
                INCLUDES="$INCLUDES -I$cf_path"
                ac_cpp="${ac_cpp} -I$cf_path"
                CFLAGS="$CFLAGS -I$cf_path"
                cf_result="yes"
                case $cf_path in
                /usr/include|/usr/include/*)
                        ;;
                *)
changequote(,)dnl
                        cf_temp=`echo $cf_path | sed -e s'%/[^/]*$%%'`
changequote([,])dnl
                        case $cf_temp in
                        */include)
                                INCLUDES="$INCLUDES -I$cf_temp"
                                ac_cpp="${ac_cpp} -I$cf_temp"
                                CFLAGS="$CFLAGS -I$cf_temp"
                                ;;
                        esac
                esac
        fi
        AC_MSG_RESULT($cf_result)
done
])dnl
