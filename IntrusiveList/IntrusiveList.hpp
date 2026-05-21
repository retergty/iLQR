/**
 * @file IntrusiveList.hpp
 * @brief 侵入式双向链表：节点带 prev/next/owner，支持 O(1) 插入与删除。
 */
#pragma once
#include <cassert>

template <typename T>
class IntrusiveList;

/**
 * @brief 侵入式链表节点基类：子类 Derive 可被挂入 IntrusiveList<Derive>。
 * @tparam Derive 派生类类型（CRTP）。
 */
template <typename Derive>
class IntrusiveListNode {
 public:
  IntrusiveListNode() = default;
  virtual ~IntrusiveListNode() = default;

 private:
  friend class IntrusiveList<Derive>;

  IntrusiveListNode* prev_{nullptr};
  IntrusiveListNode* next_{nullptr};
  IntrusiveList<Derive>* owner_{nullptr};
};

/**
 * @brief 侵入式双向链表：不分配节点内存，元素类型需继承 IntrusiveListNode<T>。
 * @tparam T 元素类型（即节点类型）。
 */
template <typename T>
class IntrusiveList {
 public:
  /** @brief 构造空链表。 */
  IntrusiveList() {
    end_.owner_ = this;
    end_.next_ = &end_;
    end_.prev_ = &end_;
    head_ = &end_;
  }
  ~IntrusiveList() = default;
  IntrusiveList(const IntrusiveList&) = delete;
  IntrusiveList& operator=(const IntrusiveList&) = delete;

  class iterator {
   public:
    explicit iterator(IntrusiveListNode<T>* now = nullptr) : now_(now) {};
    operator bool() const { return now_ != nullptr; }
    T& operator*() const { return static_cast<T&>(*now_); }
    T* operator->() const { return static_cast<T*>(now_); }

    iterator& operator++() {
      now_ = now_->next_;
      return *this;
    }

    iterator operator++(int) {
      iterator tmp = *this;
      now_ = now_->next_;
      return tmp;
    }
    iterator& operator--() {
      now_ = now_->prev_;
      return *this;
    }
    iterator operator--(int) {
      iterator tmp = *this;
      now_ = now_->prev_;
      return tmp;
    }
    bool operator==(const iterator& rhs) const { return now_ == rhs.now_; }
    bool operator!=(const iterator& rhs) const { return now_ != rhs.now_; }
    iterator operator+(int index) {
      iterator tmp = *this;
      tmp += index;
      return tmp;
    }
    iterator& operator+=(int index) {
      for (int i = 0; i < index; ++i) {
        now_ = now_->next_;
      }
      return *this;
    }

    friend class IntrusiveList<T>;

   private:
    IntrusiveListNode<T>* now_;
  };

  class const_iterator {
   public:
    explicit const_iterator(const IntrusiveListNode<T>* now = nullptr)
        : now_(now) {};
    const_iterator(const iterator& other) : now_(other.now_) {}

    operator bool() const { return now_ != nullptr; }
    const T& operator*() const { return static_cast<const T&>(*now_); }
    const T* operator->() const { return static_cast<const T*>(now_); }

    const_iterator& operator++() {
      now_ = now_->next_;
      return *this;
    }

    const_iterator operator++(int) {
      const_iterator tmp = *this;
      now_ = now_->next_;
      return tmp;
    }
    const_iterator& operator--() {
      now_ = now_->prev_;
      return *this;
    }
    const_iterator operator--(int) {
      const_iterator tmp = *this;
      now_ = now_->prev_;
      return tmp;
    }
    bool operator==(const const_iterator& rhs) const {
      return now_ == rhs.now_;
    }
    bool operator!=(const const_iterator& rhs) const {
      return now_ != rhs.now_;
    }
    const_iterator operator+(int index) {
      const_iterator tmp = *this;
      tmp += index;
      return tmp;
    }
    const_iterator& operator+=(int index) {
      for (int i = 0; i < index; ++i) {
        now_ = now_->next_;
      }
      return *this;
    }

    friend class IntrusiveList<T>;

   private:
    const IntrusiveListNode<T>* now_;
  };

  /**
   * Remove the node from List
   * @param [in] Node: point to node
   * @return -1: error, 0 success
   */
  int remove(IntrusiveListNode<T>* node_ptr) {
    assert(node_ptr != nullptr);

    IntrusiveList* owner = node_ptr->owner_;
    if (owner == nullptr) {
      return -1;
    }
    assert(node_ptr != &owner->end_);

    IntrusiveListNode<T>* prev = node_ptr->prev_;
    IntrusiveListNode<T>* next = node_ptr->next_;

    prev->next_ = next;
    next->prev_ = prev;

    // if it is head, change it
    if (owner->head_ == node_ptr) {
      owner->head_ = next;
    }

    node_ptr->owner_ = nullptr;
    node_ptr->prev_ = nullptr;
    node_ptr->next_ = nullptr;

    return 0;
  }

  /**
   * insert the node to this List ahead pos
   * @param [in] pos: insert position
   * @param [in] Node: node to be inserted
   * @return iterator point to inserted node
   */
  iterator insert(iterator pos, IntrusiveListNode<T>& node) {
    assert(pos.now_ != nullptr);
    assert(pos.now_->owner_ == this);

    if (node.owner_) {
      // this node is in other list, remove it first
      remove(&node);
    }

    IntrusiveListNode<T>* insert_ptr = pos.now_;

    node.next_ = insert_ptr;
    node.prev_ = insert_ptr->prev_;
    node.owner_ = this;

    insert_ptr->prev_->next_ = &node;
    insert_ptr->prev_ = &node;

    if (insert_ptr == head_) {
      head_ = &node;
    }

    iterator ret{&node};
    return ret;
  }

  /**
   * Erases the specified elements from the container.
   * @param [in] pos: removes the element at pos
   * @return Iterator following the last removed element.If pos refers to the
   * last element, then the end() iterator is returned.
   */
  iterator erase(iterator pos) {
    assert(pos.now_ != nullptr);
    assert(pos.now_->owner_ == this);
    assert(pos.now_ != &end_);
    IntrusiveListNode<T>* node_ptr = pos.now_;

    IntrusiveListNode<T>* last = node_ptr->prev_;
    IntrusiveListNode<T>* next = node_ptr->next_;

    last->next_ = next;
    next->prev_ = last;

    // if it is head, change it
    if (head_ == node_ptr) {
      head_ = next;
    }

    node_ptr->owner_ = nullptr;
    node_ptr->prev_ = nullptr;
    node_ptr->next_ = nullptr;

    return iterator(next);
  }

  iterator begin() { return iterator(head_); }

  iterator end() { return iterator(&end_); }

  const_iterator begin() const { return const_iterator(head_); }

  const_iterator end() const { return const_iterator(&end_); }

  const_iterator cbegin() const { return begin(); }

  const_iterator cend() const { return end(); }

 private:
  IntrusiveListNode<T>* head_;
  IntrusiveListNode<T> end_;
};