#include <Python.h>
#include <cassert>
#include <vector>
#include <map>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __i386
__inline__ uint64_t rdtsc() {
    uint64_t x;
    __asm__ volatile ("rdtsc" : "=A" (x));
    return x;
}
#elif __amd64
__inline__ uint64_t rdtsc() {
    uint64_t a, d;
    __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
    return (d<<32) | a;
}
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
    //PyObject **polList;
    PyObject **wtList;
    PyObject **cfList;
    int32_t *rem;
    int32_t *rem_weight;
    PyObject *PyZero = NULL;
    PyObject *ret;
    PyObject *p_coeffs = NULL;
    PyObject *p_weights = NULL;
    PyObject *r_coeffs = NULL;
    PyObject *r_weights = NULL;
    uint32_t Ngterm,Nbase,Nvar,Nweight,i,j,term_index;
    int32_t *htIndices;
    int32_t *p_htIndices;  /* working space */
    uint64_t t0,t1,t2,t3,t_red=0,t_tot;
    PyObject **flatten_coeffs;
    int32_t *flatten_degrees;
    int32_t *flatten_weights;
    uint32_t *startIndex;

    t0 = rdtsc();
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

    //initialize
    htList = (PyObject**)malloc(sizeof(PyObject*)*Nbase);
    wtList = (PyObject**)malloc(sizeof(PyObject*)*Nbase);
    cfList = (PyObject**)malloc(sizeof(PyObject*)*Nbase);
    PyZero = PyInt_FromLong(0);
    Py_INCREF(PyZero);
    rem = (int32_t*)malloc(sizeof(int32_t)*Nvar);
    rem_weight = (int32_t*)malloc(sizeof(int32_t)*Nweight);
    htIndices = (int32_t*)malloc(sizeof(int32_t)*Nvar*Nbase);
    p_htIndices = (int32_t*)malloc(sizeof(int32_t)*Nvar);

    p_coeffs = PyObject_GetAttrString( p , "coeffs");
    p_weights = PyObject_GetAttrString( p , "weights" );
    ret = NULL;
    Ngterm = 0;
    for(i = 0 ; i < Nbase ; i++){
       PyObject *q = PyList_GetItem(arg_polList , i);
       wtList[i] = PyObject_GetAttrString( q , "weights" );
       cfList[i] = PyObject_GetAttrString( q , "coeffs" );
       Ngterm += (uint32_t)PyDict_Size( wtList[i] );
    }
    flatten_weights = (int32_t*)malloc(sizeof(int32_t)*Ngterm*Nweight);
    flatten_degrees = (int32_t*)malloc(sizeof(int32_t)*Ngterm*Nvar);
    flatten_coeffs = (PyObject**)malloc(sizeof(PyObject*)*Ngterm);
    startIndex = (uint32_t*)malloc(sizeof(uint32_t)*Nbase);
    for(i = 0,term_index=0 ; i < Nbase ; i++){
         PyObject *q_wt = wtList[i];
         PyObject *q_cf = cfList[i];
         PyObject *q_weights_keys = PyDict_Keys( q_wt );
         uint32_t k,kmax;
         startIndex[i] = term_index;
         kmax = (uint32_t)PyDict_Size( q_wt );
         for(k = 0 ; k < kmax ; k++,term_index++){
              PyObject *q_key_k = PyList_GetItem(q_weights_keys , k);
              PyObject *q_weight_val_k = PyDict_GetItem(q_wt , q_key_k );
              PyObject *q_coeff_val_k = PyDict_GetItem(q_cf , q_key_k );
              flatten_coeffs[term_index] = q_coeff_val_k;
              for(j = 0 ; j < Nweight ; j++){
                  flatten_weights[term_index*Nweight+j] = PyInt_AsLong( PyTuple_GetItem(q_weight_val_k , j) );
              }
              for(j = 0 ; j < Nvar ; j++){
                  flatten_degrees[term_index*Nvar+j] = PyInt_AsLong( PyTuple_GetItem(q_key_k , j) );
              } 
         }
    }

    /** compute leading-term list **/
    for(i = 0 ; i < Nbase ; i++){
       uint32_t k,kmax;
       PyObject *ht;
       PyObject *hw = NULL;
       PyObject *q_weights = wtList[i];
       kmax = (uint32_t)PyDict_Size( q_weights );
       if(kmax==0)goto cleanup;
       {
           PyObject *q_weights_keys = PyDict_Keys( q_weights );
           ht = PyList_GetItem(q_weights_keys , 0);
           hw = PyDict_GetItem(q_weights , ht);
           uint32_t k_best = 0;
           for(k = 1 ; k < kmax ; k++){
              PyObject *q_weights_key_k = PyList_GetItem(q_weights_keys , k);
              PyObject *q_weights_val_k = PyDict_GetItem(q_weights , q_weights_key_k );
              if( PyObject_Compare( q_weights_val_k , hw )==1 ){
                 hw = q_weights_val_k;
                 ht = q_weights_key_k;
                 k_best = k;
              }
           }
           htList[i] = ht;
           for(j = 0 ; j < Nvar ; j++){
               int32_t q_ht_cdeg = PyInt_AsLong( PyTuple_GetItem( ht , j) );
               htIndices[(i*Nvar)+j] = q_ht_cdeg;
           }
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
        PyObject *p_ht = NULL;
        PyObject *p_hw = NULL;
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
        for(j = 0 ; j < Nvar ; j++)p_htIndices[j] = PyInt_AsLong( PyTuple_GetItem( p_ht , j) );
        Py_DECREF( p_weights_keys );

        for(i = 0; i < Nbase ; i++){
            redble = true;
            /** check if q.tip divides p.tip **/
            for(j = 0 ; j < Nvar ; j++){
               int32_t p_ht_cdeg = p_htIndices[j];
               int32_t q_ht_cdeg = htIndices[(i*Nvar)+j];
               int32_t rem_cdeg = p_ht_cdeg - q_ht_cdeg;
               if(rem_cdeg<0){
                  redble = false;
                  break;
               }
               rem[j] = rem_cdeg;
            }

            /** reduce p if there is q which tip divides p.tip **/
            if(redble){
               PyObject *q_ht = htList[i];
               PyObject *q_weights = wtList[i];
               PyObject *q_coeffs = cfList[i];
               PyObject *q_hw = PyDict_GetItem( q_weights , q_ht );
               PyObject *p_hc;
               PyObject *q_hc;
               PyObject *ratio_hc;
               uint32_t term_start,term_end;
               //t2 = rdtsc();
               for(j = 0 ; j < Nweight ; j++){
                  int32_t p_ht_cweight = PyInt_AsLong( PyTuple_GetItem( p_hw , j) );
                  int32_t q_ht_cweight = PyInt_AsLong( PyTuple_GetItem( q_hw , j) );
                  rem_weight[j] = p_ht_cweight-q_ht_cweight;
               }
#if 1
               p_hc = PyDict_GetItem( p_coeffs , p_ht );
               q_hc = PyDict_GetItem( q_coeffs , q_ht );
               ratio_hc = PyNumber_Negative( PyNumber_Divide(p_hc , q_hc) );
#else
               p_hc = PyNumber_Negative( PyDict_GetItem( p_coeffs , p_ht ) );
               q_hc = PyDict_GetItem( q_coeffs , q_ht );
#endif

               t2 = rdtsc();
               term_start = startIndex[i];
               if(i==Nbase-1){
                  term_end = Ngterm;
               }else{
                  term_end = startIndex[i+1];
               }
               for(k = term_start ; k < term_end ; k++){
                  PyObject *q_coeff_val = flatten_coeffs[k];
                  PyObject *new_key = PyTuple_New( Nvar );
                  for(j = 0 ; j < Nvar ; j++){
                     int32_t rem_cdeg = rem[j];
                     int32_t q_cdeg = flatten_degrees[k*Nvar+j];
                     PyObject *tmp = PyInt_FromLong(rem_cdeg+q_cdeg);
                     PyTuple_SetItem( new_key , j , tmp );
                  }

                  if( PyDict_Contains(p_coeffs , new_key)==1 ){
                      PyObject *old_coeff = PyDict_GetItem( p_coeffs , new_key );
#if 1
                      PyObject *tmp = PyNumber_Multiply(q_coeff_val , ratio_hc);
                      PyObject * new_coeff = PyNumber_Add( old_coeff , tmp );
#else
                      PyObject *tmp = PyNumber_Multiply(q_coeff_val , p_hc);
                      PyObject *tmp2 = PyNumber_Multiply( old_coeff , q_hc);
                      PyObject * new_coeff = PyNumber_Add( tmp2 , tmp );
                      Py_DECREF(tmp2);
#endif
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
                      //PyObject *q_weight_val = PyDict_GetItem( q_weights , q_ckey );
                      //int32_t *q_weight_val = &flatten_weights[(term_start+k)*Nweight];
#if 1
                      PyObject *new_coeff = PyNumber_Multiply(q_coeff_val , ratio_hc);
#else
                      PyObject *new_coeff = PyNumber_Multiply(q_coeff_val , p_hc);
#endif
                      PyDict_SetItem( p_coeffs , new_key , new_coeff );
                      Py_DECREF( new_coeff );
                      for(j = 0 ; j < Nweight ; j++){
                          int32_t w1 = flatten_weights[k*Nweight+j];
                          int32_t w2 = rem_weight[j];
                          PyTuple_SetItem( new_weight , j , PyInt_FromLong( w1+w2 ) );
                      }

                      PyDict_SetItem( p_weights , new_key , new_weight );
                      Py_DECREF( new_weight );
                  }
                  Py_DECREF( new_key );
               }
#if 1
               Py_DECREF( ratio_hc );
#endif
               t3 = rdtsc();
               t_red += (t3-t2);
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
           p_ht = NULL;
        }
    }
    ret = PyTuple_Pack(2 , r_coeffs , r_weights);

cleanup:
    for(i = 0 ; i < Nbase ; i++){
        if(htList[i])Py_DECREF( htList[i] );
        //Py_DECREF(wtList[i]);
        //Py_DECREF(cfList[i]);
    }
    free(wtList);
    free(cfList);
    free( htList );
    free( startIndex );
    free( flatten_degrees );
    free( flatten_weights );
    free( flatten_coeffs );
    //free( polList );
    if(htIndices!=NULL)free(htIndices);
    if(p_htIndices!=NULL)free(p_htIndices);
    if(rem!=NULL)free(rem);
    if(rem_weight!=NULL)free(rem_weight);
    if(p_coeffs!=NULL)Py_DECREF( p_coeffs );
    if(p_weights!=NULL)Py_DECREF( p_weights );
    if(r_coeffs!=NULL)Py_DECREF( r_coeffs );
    if(r_weights!=NULL)Py_DECREF( r_weights );
    if(PyZero!=NULL)Py_DECREF( PyZero );
    t1 = rdtsc();
    t_tot = (t1-t0);
    //printf("t_tot=%ld , t_red=%ld\n",t_tot,t_red);
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
