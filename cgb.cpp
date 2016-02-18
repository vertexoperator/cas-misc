#include <Python.h>
#include <cassert>

#ifdef __cplusplus
extern "C" {
#endif

inline int32_t _cx_ideg(PyObject *n , int32_t Nvar){
    int32_t ret = 0;
    PyObject *p_Nvar = PyLong_FromLong(Nvar);
    PyObject *n0 = PyLong_FromLong(1);
    PyObject *tmp;
    /* n_top = (Nvar-1)*n + 1 */
    PyObject *n_top = PyLong_FromLong(Nvar-1);
    tmp = n_top;
    n_top = PyNumber_InPlaceMultiply(tmp , n);
    Py_DECREF(tmp);
    tmp = n_top;
    n_top = PyNumber_InPlaceAdd(tmp , n0);
    Py_DECREF(tmp);
    while(1){
       if(PyObject_Compare(n0,n_top)==1)break;
       ret++;
       tmp = n0;
       n0 = PyNumber_InPlaceMultiply(tmp , p_Nvar);
       Py_DECREF(tmp);
    }
    Py_DECREF( p_Nvar );
    Py_DECREF( n0 );
    Py_DECREF( n_top );
    return (ret-1);
}


static PyObject *cx_im_mul(PyObject *self , PyObject *args){
    PyObject *n,*m,*p_Nvar;
    int32_t Nvar;
    int32_t n_deg,m_deg;
    uint32_t i;
    PyObject *buf[12];
    PyObject *p1,*p2,*base,*base1,*base2,*ret;
    if(!PyArg_UnpackTuple(args, "ref", 1, 3, &n , &m , &p_Nvar)){
        return NULL;
    }
    Nvar = PyLong_AsLong(p_Nvar);
    if(Nvar==0){
       return NULL;
    }else if(Nvar==1){
       return PyNumber_Add(n , m);
    }
    n_deg = _cx_ideg(n , Nvar);
    m_deg = _cx_ideg(m , Nvar);
    buf[0] = PyLong_FromLong(1);
    buf[1] = PyLong_FromLong(Nvar-1);
    buf[2] = PyLong_FromLong(n_deg);
    buf[3] = PyLong_FromLong(m_deg);
    p1 = PyNumber_Power(p_Nvar , buf[2] , Py_None);
    p2 = PyNumber_Power(p_Nvar , buf[3] , Py_None);
    /* base = (p1*p2-1)//(Nvar-1) */
    buf[4] = PyNumber_Multiply(p1,p2);
    buf[5] = PyNumber_Subtract(buf[4] , buf[0]);
    base = PyNumber_Divide(buf[5] , buf[1]);
    buf[6] = PyNumber_Subtract(p1 , buf[0]);
    buf[7] = PyNumber_Subtract(p2 , buf[0]);
    base1 = PyNumber_Divide( buf[6] , buf[1] );
    base2 = PyNumber_Divide( buf[7] , buf[1] );
    buf[8] = PyNumber_Subtract(m , base2);
    buf[9] = PyNumber_Subtract(n , base1);
    buf[10] = PyNumber_Multiply(buf[9] , p2);  /* p2*(n-base1) */
    buf[11] = PyNumber_Add(buf[8] , buf[10]);  /* p2*(n-base1) + (m-base2) */
    ret = PyNumber_Add(base , buf[11]);
    Py_DECREF(p1);
    Py_DECREF(p2);
    Py_DECREF(base);
    Py_DECREF(base1);
    Py_DECREF(base2);
    for(i = 0; i < 12 ; i++){
         Py_DECREF(buf[i]);
    }
    return ret;
}


static PyObject *cx_dp_nf(PyObject *self , PyObject *args){
    PyObject *p;
    PyObject *arg_polList;
    PyObject **htList;
    PyObject **polList;
    int32_t *rem;
    int32_t *rem_weight;
    PyObject *PyZero = NULL;
    PyObject *ret;
    PyObject *p_coeffs;
    PyObject *p_weights;
    PyObject *r_coeffs;
    PyObject *r_weights;
    uint32_t Nbase,Nvar,Nweight,i,j;
    
    if(!PyArg_UnpackTuple(args, "ref", 1, 2, &p , &arg_polList)){
        return NULL;
    }
    if(p==NULL || arg_polList==NULL){
        return NULL;
    }
    Nbase = PyObject_Length( arg_polList );
    if(Nbase==0)return p;
    {
        PyObject *p_Nvar,*p_Nweight;
        p_Nvar = PyObject_GetAttrString( p , "Nvar");
        p_Nweight = PyObject_GetAttrString( p , "Nweight");
        Nvar = PyInt_AsLong( p_Nvar );
        Nweight = PyInt_AsLong( p_Nweight );
        Py_DECREF( p_Nvar );
        Py_DECREF( p_Nweight );
    }

    polList = (PyObject**)malloc(sizeof(PyObject*)*Nbase);
    htList = (PyObject**)malloc(sizeof(PyObject*)*Nbase);
    PyZero = PyInt_FromLong(0);
    Py_INCREF(PyZero);
    rem = (int32_t*)malloc(sizeof(int32_t)*Nvar);
    rem_weight = (int32_t*)malloc(sizeof(int32_t)*Nweight);
    p_coeffs = PyObject_GetAttrString( p , "coeffs");
    p_weights = PyObject_GetAttrString( p , "weights" );
    ret = NULL;
    /** compute leading-term list **/
    for(i = 0 ; i < Nbase ; i++){
       uint32_t k,kmax;
       PyObject *ht;
       PyObject *max_weight = NULL;
       PyObject *q_weights = NULL;
       polList[i] = PyList_GetItem(arg_polList , i);
       q_weights = PyObject_GetAttrString(polList[i] , "weights");
       kmax = (uint32_t)PyDict_Size( q_weights );
       if( kmax==0 ){
           htList[i] = NULL;
       }else{
           PyObject *q_weights_keys = PyDict_Keys( q_weights );
           ht = PyList_GetItem(q_weights_keys , 0);
           max_weight = PyDict_GetItem(q_weights , ht);
           for(k = 1 ; k < kmax ; k++){
              PyObject *q_weights_key_k = PyList_GetItem(q_weights_keys , k);
              PyObject *q_weights_val_k = PyDict_GetItem(q_weights , q_weights_key_k );
              if( PyObject_Compare( q_weights_val_k , max_weight )==1 ){
                 max_weight = q_weights_val_k;
                 ht = q_weights_key_k;
              }
           }
           htList[i] = ht;
           Py_DECREF( q_weights_keys );
       }
       if(q_weights)Py_DECREF( q_weights );
    }

    /** normal form computation **/
    r_coeffs = PyDict_New();
    r_weights = PyDict_New();
    if(r_coeffs==NULL || r_weights==NULL)goto cleanup;
    while(1){
        uint32_t k,kmax;
        bool redble;
        /** get current leading term of p **/
        PyObject *p_ht;
        PyObject *p_hw;
        PyObject *p_weights_keys;
        kmax = (uint32_t)PyDict_Size( p_weights );
        if(kmax==0)break;
        p_weights_keys = PyDict_Keys( p_weights );
        p_ht = PyList_GetItem(p_weights_keys , 0);
        p_hw = PyDict_GetItem(p_weights , p_ht);
        for(k = 1 ; k < kmax ; k++){
           PyObject *p_weights_key_k = PyList_GetItem(p_weights_keys , k);
           PyObject *p_weights_val_k = PyDict_GetItem(p_weights , p_weights_key_k);
           if( PyObject_Compare( p_weights_val_k , p_hw )==1 ){
                p_hw = p_weights_val_k;
                p_ht = p_weights_key_k;
           }

        }

        /** reduce p if there is q which tip divides p.tip **/
        for(i = 0; i < Nbase ; i++){
            redble = true;
            PyObject *q_ht = htList[i];
            if(q_ht==NULL)continue;
            /** check if q.tip divides p.tip **/
            for(j = 0 ; j < Nvar ; j++){
               int32_t p_ht_cdeg = PyInt_AsLong( PyTuple_GetItem( p_ht , j) );
               int32_t q_ht_cdeg = PyInt_AsLong( PyTuple_GetItem( q_ht , j) );
               int32_t rem_cdeg = p_ht_cdeg - q_ht_cdeg;
               if(rem_cdeg<0){
                  redble = false;
                  break;
               }
               rem[j] = rem_cdeg;
            }

            if(redble){
               PyObject *q = polList[i];
               PyObject *q_weights_keys;
               PyObject *q_weights = PyObject_GetAttrString( q , "weights" );
               PyObject *q_coeffs = PyObject_GetAttrString( q , "coeffs" );
               PyObject *q_hw = PyDict_GetItem( q_weights , q_ht );
               PyObject *p_hc;
               PyObject *q_hc;
               PyObject *ratio_hc;

               for(j = 0 ; j < Nweight ; j++){
                  int32_t p_ht_cweight = PyInt_AsLong( PyTuple_GetItem( p_hw , j) );
                  int32_t q_ht_cweight = PyInt_AsLong( PyTuple_GetItem( q_hw , j) );
                  rem_weight[j] = p_ht_cweight-q_ht_cweight;
               }
               p_hc = PyDict_GetItem( p_coeffs , p_ht );
               q_hc = PyDict_GetItem( q_coeffs , q_ht );
               ratio_hc = PyNumber_Divide(p_hc , q_hc);
               q_weights_keys = PyDict_Keys( q_weights );
               kmax = (uint32_t)PyDict_Size( q_weights );

               for(k = 0 ; k < kmax ; k++){
                  PyObject *q_ckey = PyList_GetItem( q_weights_keys , k );
                  PyObject *q_coeff_val = PyDict_GetItem( q_coeffs , q_ckey );
                  PyObject *new_key = PyTuple_New( Nvar );
                  for(j = 0 ; j < Nvar ; j++){
                     int32_t rem_cdeg = rem[j];
                     int32_t q_cdeg = PyInt_AsLong( PyTuple_GetItem( q_ckey , j) );
                     PyObject *tmp = PyInt_FromLong(rem_cdeg+q_cdeg);
                     PyTuple_SetItem( new_key , j , tmp );
                  }

                  if( PyDict_Contains(p_coeffs , new_key)==1 ){
                      PyObject *old_coeff = PyDict_GetItem( p_coeffs , new_key );
                      PyObject *tmp = PyNumber_Multiply(q_coeff_val , ratio_hc);
                      PyObject * new_coeff = PyNumber_Subtract( old_coeff , tmp );
                      if( PyObject_Compare(new_coeff , PyZero)==0 ){
                         PyDict_DelItem( p_coeffs , new_key);
                         PyDict_DelItem( p_weights , new_key);
                         Py_DECREF( new_coeff );
                      }else{
                         PyDict_SetItem( p_coeffs , new_key , new_coeff );
                         Py_DECREF( new_coeff );
                      }
                      Py_DECREF( tmp );
                  }else{
                      PyObject *new_weight = PyTuple_New( Nweight );
                      PyObject *q_weight_val = PyDict_GetItem( q_weights , q_ckey );
                      PyObject *tmp = PyNumber_Multiply(q_coeff_val , ratio_hc);
                      PyObject *new_coeff = PyNumber_Negative( tmp );
                      PyDict_SetItem( p_coeffs , new_key , new_coeff );
                      Py_DECREF( tmp );
                      Py_DECREF( new_coeff );
                      for(j = 0 ; j < Nweight ; j++){
                          int32_t w1 = PyInt_AsLong( PyTuple_GetItem( q_weight_val , j) );
                          int32_t w2 = rem_weight[j];
                          tmp = PyInt_FromLong( w1+w2 );
                          PyTuple_SetItem( new_weight , j , tmp );
                      }

                      PyDict_SetItem( p_weights , new_key , new_weight );
                      Py_DECREF( new_weight );
                  }
                  Py_DECREF( new_key );
               }
               Py_DECREF( ratio_hc );
               Py_DECREF( q_weights_keys );
               Py_DECREF( q_weights );
               Py_DECREF( q_coeffs );

               break;
            }
        }
        /** ------- **/
        if(!redble){
           PyObject *p_hc = PyDict_GetItem( p_coeffs, p_ht );

           PyDict_SetItem( r_coeffs , p_ht , p_hc );
           PyDict_SetItem( r_weights , p_ht , p_hw );
           PyDict_DelItem( p_coeffs , p_ht);
           PyDict_DelItem( p_weights , p_ht);
        }

        Py_DECREF( p_weights_keys );
    }
    ret = PyTuple_Pack(2 , r_coeffs , r_weights);

cleanup:
    //for(i = 0 ; i < Nbase ; i++){
        //if(htList[i])Py_DECREF( htList[i] );
    //}
    free( htList );
    free( polList );
    if(rem)free(rem);
    if(rem_weight)free(rem_weight);
    Py_DECREF( p_coeffs );
    Py_DECREF( p_weights );
    Py_DECREF( r_coeffs );
    Py_DECREF( r_weights );
    Py_DECREF( PyZero );

    return ret;
}


static PyMethodDef methods[] = {
    {"cx_dp_nf" , cx_dp_nf , METH_VARARGS, ""},
    {"cx_im_mul" , cx_im_mul , METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

void initcgb(void)
{
    Py_InitModule("cgb", methods);
}

#ifdef __cplusplus
}
#endif

